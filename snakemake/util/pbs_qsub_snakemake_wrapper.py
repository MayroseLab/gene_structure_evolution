#!python
import os
import sys
from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

if "threads" in job_properties:
  ppn = job_properties['threads']
else:
  ppn = 1
if "mem_mb" in job_properties["resources"]:
  mem = ',mem=%sm' % job_properties["resources"]['mem_mb']
else:
  mem = ',mem=1g'
# job name
rule_name = job_properties["rule"]
pref = []
for x in ['species']:
  if x in job_properties["wildcards"]:
    pref.append(job_properties["wildcards"][x])
if not pref:
  base = "all_species"
else:
  base = '_'.join(pref)
queue = job_properties['params']["queue"]
priority = job_properties['params']["priority"]
stdout = job_properties['log'][0].replace('.log', '.out')
stderr = job_properties['log'][0].replace('.log', '.err')
os.system("qsub -N {base}_{rule} -p {priority} -q {queue} -l nodes=1:ppn={ppn}{mem} -o {stdout} -e {stderr} {script}".format(base=base, rule=rule_name, priority=priority, queue=queue, ppn=ppn, mem=mem, script=jobscript, stdout=stdout, stderr=stderr))
