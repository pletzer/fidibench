# LoadLeveler job
# @ shell = /bin/bash
#
# @ job_name = laplacian 
#
# @ job_type = parallel
#
# @ wall_clock_limit     = 1:00:00
#
# @ account_no = nesi00213
#
# @ network.MPI = sn_all,not_shared,US
# @ task_affinity = core(1)
#
# @ output               = $(job_name).o
# @ error                = $(job_name).e
# @ notification         = never
# @ class                = General
# @ node = 1
# @ tasks_per_node = 8
#
# @ queue

export TAU_TRACE=1
poe ./laplacian -numDims 3 -numCells 128

