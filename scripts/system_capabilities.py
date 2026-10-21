"""Read Linux CPU affinity and memory availability, including cgroup v2 limits."""
import json
import os
from pathlib import Path

GIB=1024**3

def recommend(cpu_count, affinity, total, available, quotas=(), memory_limits=()):
    cpus=min(cpu_count,affinity,*[max(1,int(q)) for q in quotas])
    capacity=total;free=available
    for limit,used in memory_limits:
        capacity=min(capacity,limit);free=min(free,max(0,limit-used))
    free=min(free,capacity)
    return {'logical_cpus':cpu_count,'available_cpus':cpus,
            'memory_total_gib':round(total/GIB,1),'memory_limit_gib':int(capacity/GIB),
            'memory_available_gib':round(free/GIB,1),'max_memory_gib':int(free/GIB),
            'recommended_threads':cpus,'recommended_memory_gib':int(free*.9/GIB),
            'memory_reserve_percent':10,'supported':free*.9>=4*GIB}

def capabilities():
    memory={}
    for line in Path('/proc/meminfo').read_text().splitlines():
        key,value=line.split(':',1);memory[key]=int(value.split()[0])*1024
    total=memory['MemTotal'];available=memory.get('MemAvailable',memory['MemFree'])
    quotas=[];limits=[];base=Path('/sys/fs/cgroup');folders={base}
    for line in Path('/proc/self/cgroup').read_text().splitlines():
        if line.startswith('0::'):
            # Namespace-relative root and every accessible parent both matter.
            path=base/line[3:].lstrip('/')
            if path.is_relative_to(base):
                while path!=base:folders.add(path);path=path.parent
    for folder in folders:
        try:
            quota,period=(folder/'cpu.max').read_text().split()
            if quota!='max':quotas.append(int(quota)/int(period))
        except (OSError,ValueError):pass
        try:
            value=(folder/'memory.max').read_text().strip()
            if value!='max':limits.append((int(value),int((folder/'memory.current').read_text())))
        except (OSError,ValueError):pass
    try:affinity=len(os.sched_getaffinity(0))
    except AttributeError:affinity=os.cpu_count() or 1
    return recommend(os.cpu_count() or 1,affinity,total,available,quotas,limits)
if __name__=='__main__':print(json.dumps(capabilities()))
