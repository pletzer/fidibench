import re
import sys
import pandas as pd
import os
import glob
import datetime
import matplotlib.pylab as plt
import seaborn as sn

plt.xticks(rotation=45)

case = "julia_vs_cxx"

#rundir = os.environ['CYLC_WORKFLOW_RUN_DIR']
rundir = f'/home/pletzera/cylc-run/{case}'
status_files = glob.glob(rundir + '/runN/log/job/1/upwind*_ncells*_m*/NN/job.status')
print(status_files)

init_times = []
exit_times = []
exec_times = []
ncells = []
job_ids = []
job_names = []
for sf in status_files:
  print(sf)
  init_time = float('inf')
  exit_time = float('inf')
  job_id = -1
  name = 'xxx'
  m = re.search(r'(upwind[^\_]+)\_ncells(\d+)', sf)
  if m:
    name = m.group(1)
    nc = int(m.group(2))
  for line in open(sf).readlines():
    m = re.match(r'CYLC_JOB_ID=(\d+)', line)
    if m:
      job_id = int(m.group(1))
    m = re.match(r'^CYLC_JOB_INIT_TIME=(\d+)\-(\d+)\-(\d+)T(\d+)\:(\d+)\:(\d+)Z', line)
    if m:
      print(line)
      y, m, d, H, M, S = int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4)), int(m.group(5)), int(m.group(6))
      init_time = datetime.datetime(year=y, month=m, day=d, hour=H, minute=M, second=S)
    m = re.match(r'^CYLC_JOB_EXIT_TIME=(\d+)\-(\d+)\-(\d+)T(\d+)\:(\d+)\:(\d+)Z', line)
    if m:
      y, m, d, H, M, S = int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4)), int(m.group(5)), int(m.group(6))
      exit_time = datetime.datetime(year=y, month=m, day=d, hour=H, minute=M, second=S)
  init_times.append(init_time)
  exit_times.append(exit_time)
  job_ids.append(job_id)
  job_names.append(name)
  ncells.append(nc)
  try:
    exec_times.append( (exit_time - init_time).seconds )
  except:
    exec_times.append(-1)

df = pd.DataFrame({'job_id': job_ids, 'exec_time': exec_times,
                   'job_name': job_names, 'ncells': ncells})
df.to_csv(f'{case}.csv')

print(df)

g = sn.barplot(data=df, x='ncells', y='exec_time', hue='job_name')
g.set_yscale("log")
plt.savefig(f'{case}.png', bbox_inches="tight")


