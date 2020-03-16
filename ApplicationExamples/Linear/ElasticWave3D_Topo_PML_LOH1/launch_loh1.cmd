#!/bin/bash
####
## optional: energy policy tags
#
# DO NOT USE environment = COPY_ALL
#@ job_type = parallel
#@ class = test 
# total_tasks = 758
#@ node = 15
#@ island_count=1
#@ tasks_per_node = 2
#@ minimize_time_to_solution = yes
#@ wall_clock_limit = 00:30:00
#@ energy_policy_tag = ExaHyPE_Euler_energy_tag
#@ job_name = loh1_pml_o4
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = ./
#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification=always
#@ notify_user=kduru@geophysik.uni-muenchen.de
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
module switch intel/18.0
module switch tbb/2018
module switch gcc/5

#for hybrid codes, uncomment this and use task_per_node=1 above
#export SHAREMEM=None
#export MP_SINGLE_THREAD=no
export OMP_NUM_THREADS=8
#export MP_PE_AFFINITY=yes
export MP_TASK_AFFINITY=core:8

poe ./ExaHyPE-Elastic ../ElasticWave3D_Topo_PML_LOH1.exahype