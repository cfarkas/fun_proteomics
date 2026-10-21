nextflow.enable.dsl = 2

include { DIA_PASEF } from './modules/diapasef/main'
params.dia_config = null

// Primary result: protein differential abundance (peptide presentation for HLA).
// Gene counts/comparisons are optional and never alter the protein evidence.
// --continue reuses matching completed Nextflow tasks and verified FragPipe searches.

params.config_file = null
params.validate_only = false
params.dry_run = false
params.outdir = null
params.help = false
params.interactively = false
params.ui_port = 8765
params.verbose = null
params.verbosity = null
params.config = null
params.threads = null
params.cpus = null
params.memory = null
params.continue = false
params.resume = false
params.report_only = false
params.rscript = System.getenv('PROTEOMICS_RSCRIPT') ?: 'auto'

def quote(value) { "'" + value.toString().replace("'", "'\\''") + "'" }
def absoluteInput(value, base) {
    def path = java.nio.file.Paths.get(value.toString())
    return (path.isAbsolute() ? path : base.resolve(path)).normalize().toString()
}

def boundedInteger(value, label, minimum, maximum) {
    if (!(value.toString() ==~ /[0-9]+/)) error "${label} must be an integer from ${minimum} to ${maximum}."
    def number = new BigInteger(value.toString())
    if (number < minimum || number > maximum) error "${label} must be between ${minimum} and ${maximum}."
    return number.intValue()
}

def prepareRuntime(project, runtime, requestedR, needR) {
    def settings = java.nio.file.Files.createTempFile('dda-pasef-launch-', '.json')
    def command = ['bash', "${project}/scripts/bootstrap_runtime.sh".toString(), runtime.toString(), requestedR.toString(), settings.toString(), needR.toString()]
    def process = new ProcessBuilder(command).inheritIO().start()
    def exitCode = process.waitFor()
    def metadata = exitCode == 0 ? settings.toFile().text : '{}'
    java.nio.file.Files.deleteIfExists(settings)
    if (exitCode != 0) error 'Dependency setup failed. See the setup output above, resolve the error and rerun the same command.'
    return new groovy.json.JsonSlurper().parseText(metadata)
}

def detectedResources(project, runtime) {
    def process = new ProcessBuilder(["${runtime}/env/bin/python".toString(), "${project}/scripts/system_capabilities.py".toString()]).start()
    def output = process.inputStream.text
    if (process.waitFor() != 0) error 'Could not detect CPU and memory capacity.'
    return new groovy.json.JsonSlurper().parseText(output)
}

def outputTree(cfg, root, validationOnly, dryRun) {
    def lines = ["[OUTPUT TREE] ${root}/", '  database/search.fasta       verified search database']
    if (validationOnly) lines << '  VALIDATION_COMPLETE.txt     input checks only'
    else {
        lines << '  analysis_details/<cohort>/<dda|ddaplus>/run/       FragPipe files and logs'
        lines << '  analysis_details/<cohort>/<dda|ddaplus>/tables_*/  detailed tables'
        if (!dryRun) {
            lines << '  analysis_details/final/     summary report and workbook'
            lines << '  final_results/'
            lines << '    START_HERE.html, results_index.csv   searchable result dashboard and test catalogue'
            lines << '    main_results/differential_expression/   primary protein results and supporting peptide tests'
            lines << '    matrices/<cohort>/<dda|ddaplus>/protein_counts.csv   primary count matrix'
            lines << '    MATRICES.html, MATRICES.csv   all count/abundance matrix links'
            lines << '    main_results/counts/   raw counts, normalized matrices and z-scores by cohort/search'
            lines << '    main_results/quality_control/   sample assignments, coverage and testing summaries'
            lines << '    FINAL_RESULTS_COMPLETE'
            if (cfg.optional_analyses?.gene_counts != false || cfg.optional_analyses?.peptide_isoforms || cfg.optional_analyses?.statistics?.enabled || cfg.optional_analyses?.drugs?.interaction?.enabled) lines << '    optional_analyses/       selected gene counts, peptide/isoform and interaction reports'
            lines << '  analysis_details/reporting/  model inputs, filters and gene-assignment audits'
        }
    }
    return lines.join('\n')
}

process VALIDATE_INPUTS {
    cache false
    input:
    path config_file
    val config_base
    path schema
    path code
    path raw_inputs, stageAs: 'raw_inputs/*/*'
    val output_root
    output:
    path 'validated.json', emit: plan
    path 'input_audit.json', emit: audit
    script:
    """
    python ${quote(code)}/pipeline.py validate --config ${quote(config_file)} --base ${quote(config_base)} --schema ${quote(schema)}
    """
}

process PREPARE_TOOLS {
    input:
    val plan
    path archives, stageAs: 'archives/*/*'
    path versions
    path code
    output:
    path 'toolset'
    script:
    """
    cat > plan.json <<'PLAN_JSON_EOF'
${plan}
PLAN_JSON_EOF
    python ${quote(code)}/pipeline.py tools --plan plan.json --versions ${quote(versions)}
    """
    stub:
    """
    mkdir -p toolset
    echo '{"stub":true}' > toolset/software_versions.json
    """
}

process BUILD_DATABASE {
    publishDir { "${output_root}/database" }, mode: 'copy', pattern: 'search.fasta'
    input:
    val plan
    path fasta_inputs, stageAs: 'fasta_inputs/*/*'
    path code
    val output_root
    output:
    path 'search.fasta', emit: fasta
    path 'database_provenance.json', emit: provenance
    script:
    """
    cat > plan.json <<'PLAN_JSON_EOF'
${plan}
PLAN_JSON_EOF
    python ${quote(code)}/pipeline.py database --plan plan.json
    """
}

process VALIDATION_COMPLETE {
    publishDir { "${output_root}" }, mode: 'copy'
    input:
    path tools
    path database
    val output_root
    output:
    path 'VALIDATION_COMPLETE.txt'
    script:
    """
    echo 'Input audit, registered archive checksums, software extraction and database validation completed. No spectra were searched.' > VALIDATION_COMPLETE.txt
    """
}

process CHECK_SEARCH {
    tag "${meta.cohort} ${meta.mode}"
    cache false
    input:
    tuple val(meta), val(plan)
    path tools
    path database
    path template
    path code
    val reuse_completed
    val cache_workdir
    output:
    tuple val(meta), val(plan), path('search_reuse.json')
    script:
    """
    cat > plan.json <<'PLAN_JSON_EOF'
${plan}
PLAN_JSON_EOF
    python ${quote(code)}/search_checkpoint.py --plan plan.json --cohort ${quote(meta.cohort)} --mode ${quote(meta.mode)} \
      --tools ${quote(tools)} --database ${quote(database)} --template ${quote(template)} --workdir ${quote(cache_workdir)} ${reuse_completed ? '--reuse' : ''}
    """
}

process FRAGPIPE {
    tag "${meta.cohort} ${meta.mode}"
    cpus { resources.cpus }
    memory { "${resources.memory_gib} GB" }
    maxForks 1
    publishDir { "${output_root}/analysis_details/${meta.cohort}/${meta.slug}" }, mode: 'copy'
    input:
    tuple val(meta), val(plan)
    path tools
    path database
    path template
    path code
    val resources
    val output_root
    output:
    tuple val(meta), path('run')
    script:
    """
    cat > plan.json <<'PLAN_JSON_EOF'
${plan}
PLAN_JSON_EOF
    python ${quote(code)}/pipeline.py search --plan plan.json \
      --cohort ${quote(meta.cohort)} --mode ${quote(meta.mode)} \
      --tools ${quote(tools)} --database ${quote(database)} --template ${quote(template)} \
      --cpus ${task.cpus} --memory ${resources.memory_gib} ${params.dry_run ? '--dry-run' : ''}
    """
    stub:
    """
    cat > plan.json <<'PLAN_JSON_EOF'
${plan}
PLAN_JSON_EOF
    python ${quote(code)}/pipeline.py stub-search --plan plan.json --cohort ${quote(meta.cohort)} --mode ${quote(meta.mode)}
    """
}

process COUNT_TABLES {
    tag "${meta.cohort} ${meta.mode}"
    publishDir { "${output_root}/analysis_details/${meta.cohort}/${meta.slug}" }, mode: 'copy'
    input:
    tuple val(meta), path(run)
    val plan
    path code
    val output_root
    output:
    path "tables_${meta.cohort}_${meta.slug}"
    script:
    """
    cat > plan.json <<'PLAN_JSON_EOF'
${plan}
PLAN_JSON_EOF
    python ${quote(code)}/reporting.py report --plan plan.json --run ${quote(run)} \
      --cohort ${quote(meta.cohort)} --mode ${quote(meta.mode)} --out tables_${meta.cohort}_${meta.slug}
    """
}

process FINAL_REPORT {
    publishDir { "${output_root}/analysis_details" }, mode: 'copy'
    input:
    path reports
    val plan
    path code
    val output_root
    output:
    path 'final'
    script:
    def inputs = reports.collect { reportPath -> quote(reportPath) }.join(' ')
    """
    cat > plan.json <<'PLAN_JSON_EOF'
${plan}
PLAN_JSON_EOF
    python ${quote(code)}/reporting.py finalize --plan plan.json --reports ${inputs}
    """
}

process FINAL_RESULTS {
    stageInMode 'copy'
    publishDir { "${output_root}" }, mode: 'copy'
    cpus { Math.min(resources.cpus as int, 8) }
    memory '16 GB'
    input:
    path reports
    path search_summary
    val plan
    path code
    val resources
    val output_root
    val rscript
    output:
    path 'final_results'
    path 'analysis_details/reporting'
    script:
    def inputs = reports.collect { quote(it) }.join(' ')
    def summaryOption = search_summary ? '--search-summary '+quote(search_summary instanceof List ? search_summary[0] : search_summary) : ''
    """
    cat > plan.json <<'PLAN_JSON_EOF'
${plan}
PLAN_JSON_EOF
    python ${quote(code)}/final_results.py --plan plan.json --reports ${inputs} \
      --out final_results --rscript ${quote(rscript)} --rlib ${quote(projectDir.toString()+'/.runtime/R-library')} ${summaryOption}
    """
}

workflow {
    def progressMonitor = null
    def publishedMetadata = []
    def legacyProvenanceDir = null
    workflow.onComplete = {
        // Internal metadata remains in work/ for reporting and resume, not published results.
        publishedMetadata.each { metadata ->
            java.nio.file.Files.deleteIfExists(java.nio.file.Paths.get(metadata.toString()))
        }
        if (legacyProvenanceDir != null) {
            def oldDirectory = new File(legacyProvenanceDir.toString())
            if (oldDirectory.isDirectory() && oldDirectory.list().length == 0) oldDirectory.delete()
        }
        if (progressMonitor != null) {
            progressMonitor.destroy()
            progressMonitor.waitFor(3, java.util.concurrent.TimeUnit.SECONDS)
        }
        if (params.interactively.toString() == 'true') {
            log.info "[INTERFACE ${workflow.success ? 'CLOSED' : 'FAILED'}] ${workflow.runName} | exit ${workflow.exitStatus}. Analysis status is recorded separately in the output logs; reopen the interface to reconnect."
        } else {
            log.info "[WORKFLOW ${workflow.success ? 'COMPLETE' : 'FAILED'}] ${workflow.runName} | elapsed ${workflow.duration} | exit ${workflow.exitStatus}"
            if (!workflow.success) log.info "Inspect the Nextflow controller log and the failing task's .command.err / run/pipeline.log. After resolving the error, use -resume ${workflow.runName} with the same configuration and work directory."
        }
    }
    if (params.help) {
        log.info file("${projectDir}/assets/cli_help.txt").text
        return
    }
    def allowed = ['config_file', 'config', 'validate_only', 'dry_run', 'outdir', 'help', 'interactively', 'ui_port',
                   'verbose', 'verbosity', 'threads', 'cpus', 'memory', 'runtime_dir', 'continue', 'resume', 'report_only', 'rscript', 'dia_config']
    def unknown = params.keySet().findAll { !allowed.contains(it) }
    if (unknown) error "Unknown pipeline option(s): ${unknown.collect { '--' + it }.join(', ')}. See --help."
    if (params.verbosity != null && params.verbose != null) error 'Use --verbosity OR legacy --verbose, not both.'
    if (params.verbose != null && !['true', 'false'].contains(params.verbose.toString())) error '--verbose accepts true or false; prefer --verbosity 1–4.'
    def verbosity = params.verbosity != null ? boundedInteger(params.verbosity, '--verbosity', 1, 4) : params.verbose != null ? (params.verbose.toString() == 'true' ? 4 : 1) : 3
    if (params.threads != null && params.cpus != null) error 'Use --threads OR --cpus, not both.'
    if (params.config != null && params.config_file != null) error 'Use --config OR --config_file, not both.'
    // A folder selects portable continuation; bare --continue keeps the DDA resume interface.
    if (params.continue != null && !['true', 'false'].contains(params.continue.toString())) {
        if (params.config != null || params.config_file != null || params.dia_config != null) error 'Use --continue <run-folder> alone; its saved settings identify the analysis.'
        def command = [System.getenv('PROTEOMICS_PYTHON') ?: 'python3',
            "${projectDir}/scripts/continue_folder.py".toString(), '--continue', params.continue.toString()]
        if (params.outdir != null) command += ['--outdir', params.outdir.toString()]
        if (params.rscript != 'auto') command += ['--rscript', params.rscript.toString()]
        def continuation = new ProcessBuilder(command).inheritIO().start()
        if (continuation.waitFor() != 0) error 'Folder continuation failed; see the integrity or dependency error above.'
        return
    }
    if (params.dia_config != null) {
        if (params.config != null || params.config_file != null) error 'Use --dia_config alone for DIA-PASEF, or --config for DDA-PASEF.'
        if (params.interactively || params.report_only || params.dry_run || params.continue || params.resume) error 'DIA-PASEF supports --validate_only and native -resume. Browser setup, --report_only, --dry_run and --continue are not supported by this branch.'
        if (nextflow.version.toString() != '26.04.2') error 'Set NXF_VER=26.04.2 for the tested DIA-NN workflow.'
        if (workflow.profile.split(',').contains('test') && !workflow.stubRun) error 'The test profile requires -stub-run for DIA-PASEF.'
        def diaConfig = file(params.dia_config, checkIfExists: true).toAbsolutePath()
        def command = [System.getenv('PROTEOMICS_DIA_PYTHON') ?: 'python3',
            "${projectDir}/scripts/diapasef_config.py".toString(), '--config', diaConfig.toString(),
            '--module-dir', "${projectDir}/modules/diapasef".toString()]
        if (params.threads != null || params.cpus != null) command += ['--threads', (params.threads ?: params.cpus).toString()]
        if (params.memory != null) command += ['--memory', params.memory.toString().trim().replaceFirst(/(?i)\s*(GiB|GB)$/, '').trim()]
        if (params.outdir != null) command += ['--outdir', absoluteInput(params.outdir, java.nio.file.Paths.get(launchDir.toString()))]
        def check = new ProcessBuilder(command).redirectError(ProcessBuilder.Redirect.INHERIT).start()
        def json = check.inputStream.text
        if (check.waitFor() != 0) error 'DIA-PASEF configuration validation failed; no search was started.'
        def settings = new groovy.json.JsonSlurper().parseText(json)
        settings.validate_only = params.validate_only
        settings.test_mode = workflow.stubRun && workflow.profile.split(',').contains('test')
        log.info "[DIA-PASEF] ${workflow.stubRun ? 'SYNTHETIC STUB' : params.validate_only ? 'VALIDATION ONLY' : 'DIA-NN 2.7.0'} | ${settings.expected_samples} samples | ${settings.cpus} CPUs | ${settings.memory_gib} GiB per engine task"
        settings.groups.each { type, group -> log.info "[DIA-PASEF] ${type}: ${group.samples} samples; ${group.library_references} library references" }
        log.info "[OUTPUT] ${settings.outdir} | resume=${workflow.resume}"
        DIA_PASEF(settings)
        return
    }
    def configPath = params.config ?: params.config_file ?: "${projectDir}/config.local.json"
    if (!['true', 'false'].contains(params.continue.toString())) error '--continue accepts true or false.'
    if (params.resume) error 'Use --continue, or native Nextflow -resume (one dash).'
    if (!file(configPath).exists() && params.interactively.toString() == 'true' && params.config == null && configPath.toString() == "${projectDir}/config.local.json") {
        file(configPath).text = file("${projectDir}/assets/empty_study.json").text
    }
    if (!file(configPath).exists()) error "Configuration not found: ${configPath}. Create config.local.json from examples/run.json, or pass --config_file."
    if (!['true', 'false'].contains(params.interactively.toString())) error '--interactively accepts true or false.'
    def launchRuntime = prepareRuntime(projectDir, params.runtime_dir, params.rscript, params.interactively.toString() == 'true' || (!params.validate_only && !params.dry_run && !workflow.stubRun))
    if (params.interactively.toString() == 'true') {
        if (params.validate_only || params.dry_run || workflow.stubRun) error 'Use the interface buttons to validate or run; do not combine --interactively with validation/dry-run/stub flags.'
        def port = boundedInteger(params.ui_port, '--ui_port', 1024, 65535)
        def python = new File("${params.runtime_dir}/env/bin/python")
        if (!python.isFile()) error 'Runtime setup did not produce Python.'
        def limits = detectedResources(projectDir, params.runtime_dir)
        def options = [:]
        if (params.outdir != null) options.outdir = absoluteInput(params.outdir, java.nio.file.Paths.get(launchDir.toString()))
        if (params.threads != null || params.cpus != null) options.cpus = boundedInteger(params.threads ?: params.cpus, '--threads', 1, limits.available_cpus)
        if (params.memory != null) options.memory_gib = boundedInteger(params.memory.toString().trim().replaceFirst(/(?i)\s*(GiB|GB)$/, '').trim(), '--memory', 4, limits.max_memory_gib)
        def command = [python.toString(), "${projectDir}/scripts/config_server.py".toString(), '--config', file(configPath).toAbsolutePath().toString(), '--port', port.toString(), '--detach', '--open-browser', '--overrides', groovy.json.JsonOutput.toJson(options)]
        log.info "[INTERFACE] Opening http://127.0.0.1:${port} — select inputs and conditions, then click Run analysis. No analysis starts automatically."
        log.info '[INTERFACE] This terminal displays analysis progress. Closing it leaves the setup page and analysis available. Use Stop run in the browser to stop an analysis.'
        def builder = new ProcessBuilder(command).directory(new File(projectDir.toString())).inheritIO()
        builder.environment().put('PROTEOMICS_RUNTIME_DIR', params.runtime_dir.toString())
        builder.environment().put('PROTEOMICS_RSCRIPT', launchRuntime.rscript.toString())
        builder.environment().put('PATH', "${params.runtime_dir}/env/bin:${System.getenv('PATH')}".toString())
        builder.environment().put('JAVA_HOME', "${params.runtime_dir}/env".toString())
        def server = builder.start()
        def shutdown = new Thread({ server.destroy() } as Runnable)
        Runtime.runtime.addShutdownHook(shutdown)
        def exitCode = server.waitFor()
        Runtime.runtime.removeShutdownHook(shutdown)
        if (exitCode != 0) error "Interface exited with code ${exitCode}. See the terminal output above."
        return
    }
    if (params.validate_only && params.dry_run) error 'Choose validate_only or dry_run, not both.'
    def config = file(configPath, checkIfExists: true)
    def base = config.toAbsolutePath().parent
    def cfg = new groovy.json.JsonSlurper().parseText(config.text)
    if (!cfg.samples || !cfg.database || !cfg.tools || !cfg.outdir) error 'Configuration must contain samples, database, tools and outdir.'
    def outputRoot = params.outdir ? absoluteInput(params.outdir, java.nio.file.Paths.get(launchDir.toString())) : absoluteInput(cfg.outdir, base)
    if (params.report_only) {
        cfg.outdir = outputRoot
        cfg.database.fasta = absoluteInput(cfg.database.fasta, base)
        def reports = cfg.samples.collect { it.cohort }.unique().collectMany { cohort ->
            cfg.modes.collect { mode ->
                def slug = mode == 'DDA+' ? 'ddaplus' : 'dda'
                def path = "${outputRoot}/analysis_details/${cohort}/${slug}/tables_${cohort}_${slug}"
                if (!file(path).exists()) path = "${outputRoot}/${cohort}/${slug}/tables_${cohort}_${slug}"
                file(path, checkIfExists: true)
            }
        }
        log.info '[REPORT ONLY] Rebuilding final_results from completed tables; no spectra searched.'
        def summary = file("${outputRoot}/analysis_details/final")
        FINAL_RESULTS(Channel.value(reports), Channel.value(summary.exists() ? summary : []), groovy.json.JsonOutput.toJson(cfg), file("${projectDir}/pipeline_lib"), cfg.resources ?: [cpus: 1], outputRoot, launchRuntime.rscript)
        return
    }
    legacyProvenanceDir = "${outputRoot}/provenance"
    ['validated.json', 'input_audit.json', 'VALIDATION_COMPLETE.txt'].each { filename ->
        publishedMetadata << "${legacyProvenanceDir}/${filename}"
    }
    publishedMetadata << "${outputRoot}/database/database_provenance.json"
    publishedMetadata << "${outputRoot}/analysis_details/final/provenance.json"
    cfg.samples.collect { it.cohort }.unique().each { cohort ->
        cfg.modes.each { mode ->
            def slug = mode == 'DDA+' ? 'ddaplus' : 'dda'
            // Search provenance and checkpoints stay with completed runs for --continue.
            publishedMetadata << "${outputRoot}/analysis_details/${cohort}/${slug}/tables_${cohort}_${slug}/report_provenance.json"
        }
    }
    def code = file("${projectDir}/pipeline_lib")
    def rawInputs = cfg.samples.collect { sample -> file(absoluteInput(sample.raw_path, base), checkIfExists: true) }
    def archives = cfg.tools.values().collect { archive -> file(absoluteInput(archive, base), checkIfExists: true) }
    def fastaInputs = [file(absoluteInput(cfg.database.fasta, base), checkIfExists: true)]
    if (cfg.database.contaminants_zip) fastaInputs << file(absoluteInput(cfg.database.contaminants_zip, base), checkIfExists: true)
    def limits = detectedResources(projectDir, params.runtime_dir)
    def resources = new LinkedHashMap(cfg.resources ?: [cpus: 32, memory_gib: 128])
    resources.cpus = boundedInteger(params.threads != null ? params.threads : params.cpus != null ? params.cpus : resources.cpus, '--threads / resources.cpus', 1, limits.available_cpus)
    def memoryValue = params.memory != null ? params.memory.toString().trim().replaceFirst(/(?i)\s*(GiB|GB)$/, '').trim() : resources.memory_gib
    resources.memory_gib = boundedInteger(memoryValue, '--memory / resources.memory_gib', 4, limits.max_memory_gib)
    log.info "[RUN] ${workflow.runName} | verbosity=${verbosity} | ${resources.cpus} CPUs | ${resources.memory_gib} GiB | resume=${workflow.resume}"
    log.info "[OUTPUT] ${outputRoot}"
    if (verbosity >= 2) log.info outputTree(cfg, outputRoot, params.validate_only, params.dry_run)
    if (verbosity >= 2) {
        def cohorts = cfg.samples.groupBy { sample -> sample.cohort }
        log.info """
================ DDA-PASEF execution plan ================
Run name : ${workflow.runName}
Config   : ${config.toAbsolutePath()}
Work     : ${workflow.workDir}
Output   : ${outputRoot}
Continue : ${params.continue.toString() == 'true' ? 'ON — reuse verified completed searches; run missing or changed steps' : workflow.resume ? 'ON — matching completed tasks will be CACHED' : 'OFF — use --continue to reuse completed work'}
Primary  : ${cfg.proteomics?.analysis == 'immunopeptidomics' ? 'HLA peptide analysis' : 'Protein differential abundance; gene comparisons are optional'}
Mode     : ${workflow.stubRun ? 'SYNTHETIC STUB — not scientific results' : params.validate_only ? 'VALIDATION ONLY' : params.dry_run ? 'COMMAND DRY RUN — no spectra searched' : 'FULL ANALYSIS'}
Resources: ${resources.cpus} CPUs / ${resources.memory_gib} GiB per search; searches run sequentially

[1] VALIDATE_INPUTS — configuration, samples, raw metadata and reference audit (always reruns)
[2] PREPARE_TOOLS — verify registered archives and extract pinned tools
[3] BUILD_DATABASE — verify/build target, contaminant and decoy database
"""
        if (params.validate_only) {
            log.info '[4] VALIDATION_COMPLETE — write validation marker; no searches or count tables'
        } else {
            cohorts.each { cohort, samples ->
                cfg.modes.each { mode ->
                    log.info "[4] FRAGPIPE — ${cohort} / ${mode}: ${samples.size()} samples (${verbosity >= 3 ? samples.collect { it.sample_id }.join(', ') : 'sample IDs shown at verbosity 3–4'})"
                }
            }
            log.info(cfg.proteomics?.analysis == 'immunopeptidomics' ? '    HLA search: nonspecific MSFragger → MSBooster/Percolator → 1% PSM FDR → IonQuant (MBR off); peptide counts and peptide-only differential tests' : '    Each search: independent raw staging → MSFragger → MSBooster → Percolator → protein inference/filtering → IonQuant MBR/MaxLFQ')
            if (cfg.optional_analyses?.peptide_isoforms) log.info '    Optional peptide/isoform evidence: sequence mapping against every FASTA target, with I/L equivalence'
            if (cfg.optional_analyses?.gene_counts != false && cfg.differential?.enabled && cfg.proteomics?.analysis != 'immunopeptidomics') log.info '    Optional gene counts: unambiguous assignments; separate models and FDR; protein counts unchanged'
            if (cfg.optional_analyses?.drugs?.interaction?.enabled) log.info '    Optional drug interaction: combination - drug A - drug B + control; separate from viability synergy'
            if (!params.dry_run) {
                log.info '[5] COUNT_TABLES — per-cohort/mode spectral counts, intensities and detection summaries'
                log.info '[6] FINAL_REPORT / FINAL_RESULTS — protein analyses first, supporting peptides, then selected optional gene/modification/interaction analyses'
            }
        }
        log.info 'Live task status and FragPipe logs follow. A heartbeat means no task exit code yet, not proof of scientific progress.'
        log.info 'Keep .nextflow/, the work directory and completed search outputs. --continue verifies matching searches before reusing them.'
        log.info '=============================================================='
    }
    if (verbosity >= 1) {
        def runtimePython = file("${params.runtime_dir}/env/bin/python")
        def python = runtimePython.exists() ? runtimePython.toString() : 'python3'
        progressMonitor = new ProcessBuilder(python, '-u', "${projectDir}/scripts/stream_progress.py",
            '--parent', ProcessHandle.current().pid().toString(),
            '--work-dir', workflow.workDir.toString(),
            '--since', (System.currentTimeMillis() / 1000.0).toString(),
            '--verbosity', verbosity.toString()).inheritIO().start()
    }
    VALIDATE_INPUTS(config, base.toString(), file("${projectDir}/assets/config.schema.json"), code, rawInputs, outputRoot)
    def plan = VALIDATE_INPUTS.out.plan.map { planFile ->
        def validated = new groovy.json.JsonSlurper().parseText(planFile.text)
        if (validated.resources == resources && !params.outdir) return planFile.text
        validated.resources = resources
        if (params.outdir) validated.outdir = outputRoot
        return groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(validated)) + '\n'
    }
    PREPARE_TOOLS(plan, archives, file("${projectDir}/assets/tool_versions.json"), code)
    BUILD_DATABASE(plan, fastaInputs, code, outputRoot)
    if (params.validate_only) {
        VALIDATION_COMPLETE(PREPARE_TOOLS.out, BUILD_DATABASE.out.fasta, outputRoot)
    } else {
        def jobs = plan.flatMap { planText ->
            cfg.samples.collect { sample -> sample.cohort }.unique().collectMany { cohort ->
                cfg.modes.collect { mode ->
                    tuple([cohort: cohort, mode: mode, slug: mode == 'DDA+' ? 'ddaplus' : 'dda'], planText)
                }
            }
        }
        def searchTemplate = cfg.proteomics?.analysis == 'immunopeptidomics' ? "${projectDir}/assets/fragpipe_hla.workflow" :
                             cfg.proteomics?.analysis == 'phosphoproteomics' ? "${projectDir}/modules/phosphoproteomics/fragpipe.workflow" : "${projectDir}/assets/fragpipe_lfq_mbr.workflow"
        CHECK_SEARCH(jobs, PREPARE_TOOLS.out, BUILD_DATABASE.out.fasta, file(searchTemplate), code,
                     params.continue.toString() == 'true' && !params.dry_run && !workflow.stubRun, workflow.workDir.toString())
        def decisions = CHECK_SEARCH.out.map { meta, planText, decision ->
            tuple(meta, planText, new groovy.json.JsonSlurper().parseText(decision.text).run)
        }.branch { meta, planText, cached ->
            reused: cached != null
            pending: true
        }
        def pending = decisions.pending.map { meta, planText, cached -> tuple(meta, planText) }
        def reused = decisions.reused.map { meta, planText, cached -> tuple(meta, file(cached, checkIfExists:true)) }
        FRAGPIPE(pending, PREPARE_TOOLS.out, BUILD_DATABASE.out.fasta,
                 file(searchTemplate), code, resources, outputRoot)
        def searches = FRAGPIPE.out.mix(reused)
        if (!params.dry_run) {
            COUNT_TABLES(searches, plan, code, outputRoot)
            FINAL_REPORT(COUNT_TABLES.out.collect(), plan, code, outputRoot)
            if (!workflow.stubRun) FINAL_RESULTS(COUNT_TABLES.out.collect(), FINAL_REPORT.out, plan, code, resources, outputRoot, launchRuntime.rscript)
        }
    }
}
