// Native DIA-NN tasks share main.nf's executor, task cache and resume session.
def sq(value) { "'" + value.toString().replace("'", "'\\''") + "'" }

process DIA_ENVIRONMENT {
    label 'dia_light'
    cpus { settings.test_mode ? 1 : 2 }
    memory { settings.test_mode ? '512 MB' : '4 GB' }
    input:
    val settings
    path sources
    path requirements, stageAs: 'requirements-lock.txt'
    path install_cache, stageAs: 'install_cache'
    output:
    path 'environment', emit: runtime
    script:
    """
    ${sq(settings.python)} 00_env.py --root environment --requirements requirements-lock.txt --cache-dir install_cache
    """
    stub:
    """
    mkdir -p environment
    echo 'STUB ONLY' > environment/STUB.txt
    """
}

process DIA_PREPARE_SAMPLES {
    label 'dia_light'
    cpus { settings.test_mode ? 1 : 2 }
    memory { settings.test_mode ? '512 MB' : '4 GB' }
    publishDir { settings.outdir }, mode: 'copy'
    input:
    val settings
    path runtime, stageAs: 'environment'
    path sources
    path sheet, stageAs: 'samplesheet.tsv'
    path raw_data, stageAs: 'raw_data'
    path canonical, stageAs: 'targets.fasta'
    path contaminants, stageAs: 'contaminants.fasta'
    output:
    path 'prepared/samples.tsv', emit: samples
    path 'prepared/groups/*.tsv', emit: groups
    path 'prepared/search.fasta', emit: fasta
    script:
    def contaminant_option = contaminants ? '--contaminant-fasta contaminants.fasta' : ''
    def service_types = groovy.json.JsonOutput.toJson(settings.service_matrices.keySet().toList())
    """
    environment/.venv/bin/python 01_prepare_samples.py --samplesheet samplesheet.tsv \
      --raw-root raw_data --canonical-fasta targets.fasta ${contaminant_option} \
      --expected-samples ${settings.expected_samples} --service-types ${sq(service_types)}
    """
    stub:
    """
    ${sq(settings.python)} stub_prepare.py
    """
}

process DIA_PREDICT_LIBRARY {
    label 'dia_engine'
    cpus { settings.test_mode ? 1 : settings.cpus }
    memory { settings.test_mode ? '512 MB' : "${settings.memory_gib} GB" }
    maxForks 1
    publishDir { "${settings.outdir}/analysis_details/prediction" }, mode: 'copy'
    input:
    val settings
    path runtime, stageAs: 'environment'
    path sources
    path method, stageAs: 'method.toml'
    path fasta, stageAs: 'search.fasta'
    output:
    path 'predicted', emit: prediction
    script:
    """
    environment/.venv/bin/python 02_predict_library.py --environment environment \
      --method method.toml --fasta search.fasta --threads ${task.cpus}
    """
    stub:
    """
    mkdir -p predicted
    echo 'STUB ONLY' > predicted/STUB.txt
    """
}

process DIA_BUILD_EMPIRICAL_LIBRARY {
    tag { sample_type }
    label 'dia_engine'
    cpus { settings.test_mode ? 1 : settings.cpus }
    memory { settings.test_mode ? '512 MB' : "${settings.memory_gib} GB" }
    maxForks 1
    publishDir path: { "${settings.outdir}/analysis_details/libraries/${sample_type}" }, mode: 'copy', pattern: 'library'
    input:
    val settings
    path runtime, stageAs: 'environment'
    path sources
    path method, stageAs: 'method.toml'
    path fasta, stageAs: 'search.fasta'
    path prediction, stageAs: 'predicted'
    tuple val(sample_type), path(samples, stageAs: 'samples.tsv')
    path raw_data, stageAs: 'raw_data'
    output:
    tuple val(sample_type), path('samples.tsv'), path('library'), emit: library
    script:
    """
    environment/.venv/bin/python 03_build_library.py --environment environment \
      --method method.toml --fasta search.fasta --samples samples.tsv --prediction predicted --threads ${task.cpus}
    """
    stub:
    """
    mkdir -p library
    echo 'STUB ONLY' > library/STUB.txt
    """
}

process DIA_QUANTIFY_SAMPLES {
    tag { sample_type }
    label 'dia_engine'
    cpus { settings.test_mode ? 1 : settings.cpus }
    memory { settings.test_mode ? '512 MB' : "${settings.memory_gib} GB" }
    maxForks 1
    publishDir path: { "${settings.outdir}/analysis_details/search/${sample_type}" }, mode: 'copy', pattern: 'cohort'
    input:
    val settings
    path runtime, stageAs: 'environment'
    path sources
    path method, stageAs: 'method.toml'
    path fasta, stageAs: 'search.fasta'
    tuple val(sample_type), path(samples, stageAs: 'samples.tsv'), path(library, stageAs: 'library')
    path raw_data, stageAs: 'raw_data'
    output:
    tuple val(sample_type), path('samples.tsv'), path('cohort'), emit: cohort
    script:
    """
    environment/.venv/bin/python 04_quantify_samples.py --environment environment \
      --method method.toml --fasta search.fasta --samples samples.tsv --library library --threads ${task.cpus}
    """
    stub:
    """
    mkdir -p cohort
    echo 'STUB ONLY' > cohort/STUB.txt
    """
}

process DIA_COUNT_PROTEINS {
    tag { sample_type }
    label 'dia_light'
    cpus { settings.test_mode ? 1 : 2 }
    memory { settings.test_mode ? '512 MB' : '4 GB' }
    publishDir { "${settings.outdir}/analysis_details/groups" }, mode: 'copy'
    input:
    val settings
    path runtime, stageAs: 'environment'
    path sources
    path method, stageAs: 'method.toml'
    path fasta, stageAs: 'search.fasta'
    tuple val(sample_type), path(samples, stageAs: 'samples.tsv'), path(cohort, stageAs: 'cohort'), path(service, stageAs: 'service.tsv')
    output:
    tuple val(sample_type), path("${sample_type}"), emit: results
    script:
    def service_option = service ? '--service-matrix service.tsv' : ''
    def differential = groovy.json.JsonOutput.toJson(settings.differential)
    """
    environment/.venv/bin/python 05_count_proteins.py --samples samples.tsv \
      --report cohort/report.parquet --method method.toml --outdir ${sq(sample_type)} --fasta search.fasta --differential ${sq(differential)} ${service_option}
    """
    stub:
    """
    mkdir -p ${sq(sample_type)}
    echo 'STUB ONLY - no counts' > ${sq(sample_type + '/STUB.txt')}
    """
}

process DIA_COLLECT_RESULTS {
    label 'dia_light'
    cpus { settings.test_mode ? 1 : 2 }
    memory { settings.test_mode ? '512 MB' : '4 GB' }
    publishDir { settings.outdir }, mode: 'copy'
    input:
    val settings
    path runtime, stageAs: 'environment'
    path sources
    path samples, stageAs: 'all_samples.tsv'
    path groups, stageAs: 'group_results/*'
    output:
    path 'final_results', emit: results
    script:
    def directories = [groups].flatten().collect { sq(it) }.join(' ')
    """
    environment/.venv/bin/python 06_collect_results.py --samples all_samples.tsv --groups ${directories}
    """
    stub:
    """
    mkdir -p final_results
    echo 'STUB ONLY - no scientific results' > final_results/STUB.txt
    """
}

process DIA_VALIDATION_COMPLETE {
    cpus 1
    memory '1 GB'
    publishDir { settings.outdir }, mode: 'copy'
    input:
    val settings
    path samples
    path fasta
    output:
    path 'DIA_VALIDATION_COMPLETE.txt'
    script:
    """
    echo 'DIA-PASEF samples, service matrices, FASTA and pinned DIA-NN environment validated. No spectra searched.' > DIA_VALIDATION_COMPLETE.txt
    """
}

workflow DIA_PASEF {
    take:
    settings
    main:
    def names = ['00_env.py', '01_prepare_samples.py', '02_predict_library.py', '03_build_library.py',
                 '04_quantify_samples.py', '05_count_proteins.py', '06_collect_results.py', 'common.py', 'stub_prepare.py']
    def sources = Channel.value(names.collect { file("${projectDir}/modules/diapasef/${it}", checkIfExists: true) } + ["quantification_config.py", "quantification.py", "quantification_analysis.py", "intensity_differential.R", "matrix_catalog.py"].collect { file("${projectDir}/pipeline_lib/${it}", checkIfExists: true) })
    def method = Channel.value(file(settings.method_config, checkIfExists: true))
    def raws = Channel.value(file(settings.raw_root, checkIfExists: true))
    def services = settings.service_matrices.collectEntries { key, value -> [key, file(value, checkIfExists: true)] }
    def runtime = DIA_ENVIRONMENT(settings, sources,
        Channel.value(file(settings.requirements_lock, checkIfExists: true)),
        Channel.value(file(settings.install_cache, checkIfExists: true)))
    def prepared = DIA_PREPARE_SAMPLES(settings, runtime.runtime, sources,
        Channel.value(file(settings.samplesheet, checkIfExists: true)), raws,
        Channel.value(file(settings.canonical_fasta, checkIfExists: true)),
        Channel.value(settings.contaminant_fasta ? file(settings.contaminant_fasta, checkIfExists: true) : []))
    if (settings.validate_only) {
        DIA_VALIDATION_COMPLETE(settings, prepared.samples, prepared.fasta)
    } else {
        def predicted = DIA_PREDICT_LIBRARY(settings, runtime.runtime, sources, method, prepared.fasta)
        def groups = prepared.groups.flatten().map { sheet -> tuple(sheet.baseName, sheet) }
        def libraries = DIA_BUILD_EMPIRICAL_LIBRARY(settings, runtime.runtime, sources, method, prepared.fasta, predicted.prediction, groups, raws)
        def quantified = DIA_QUANTIFY_SAMPLES(settings, runtime.runtime, sources, method, prepared.fasta, libraries.library, raws)
        def with_service = quantified.cohort.map { type, sheet, cohort -> tuple(type, sheet, cohort, services[type] ?: []) }
        def counted = DIA_COUNT_PROTEINS(settings, runtime.runtime, sources, method, prepared.fasta, with_service)
        def results = counted.results.map { type, directory -> directory }.collect()
        DIA_COLLECT_RESULTS(settings, runtime.runtime, sources, prepared.samples, results)
    }
}
