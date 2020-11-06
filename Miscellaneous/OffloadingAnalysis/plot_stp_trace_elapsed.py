#!/usr/bin/env python

import parse_stp_trace as pt
import sys
import matplotlib.pyplot as plt

#import tikzplotlib 

plt.style.use("ggplot")

#chart_name="Time Spent in STPs per Rank"
x_label="Ranks"
y_label="Total CPU Time in Space Time Predictor"

basename1 = "exahype_solvers_ADERDGSolver_PredictionJob_run_rank_"
basename2 = "exahype_solvers_ADERDGSolver_OwnMigratableJob_run_rank_"
basename3 = "exahype_solvers_ADERDGSolver_AlienMigratableJob_run_rank_"


def plot_stp_trace_single(ax,x, iterations, time, x_label, y_label, color, label_add, bottom):
  ax.set_xlabel(x_label)
  ax.set_ylabel(y_label)
  #ax.title(chart_name)
  ax.bar(x,time,0.8,label=label_add, bottom=bottom, color=color)
  #plt.bar(ranks,values_fts,0.8,Color="Red",bottom=values_stp,label="FTS")
  #plt.get_xaxis().set_ticks([0,2,4,6,8,10,12,14,16,18,20,22,24,26])
  #ax.set_ylim(0,600)

if __name__=="__main__":
  dirname = sys.argv[1]
  ranks = int(sys.argv[2])
  threads = int(sys.argv[3])
  timestep = int(sys.argv[4])

  f, ax1 = plt.subplots(1)

  (x, iterations_pred, time_pred) = pt.parseDir(dirname, ranks, threads, timestep, basename1)
  plot_stp_trace_single(ax1, x, iterations_pred, time_pred, x_label, y_label, 'b', 'Skeleton Job', 0)

  (x, iterations_own, time_own) = pt.parseDir(dirname, ranks, threads, timestep, basename2)
  plot_stp_trace_single(ax1, x, iterations_own, time_own, x_label, y_label, 'g', 'Migratable Job (local)', time_pred)

  (x, iterations_remote, time_remote) = pt.parseDir(dirname, ranks, threads, timestep, basename3)
  plot_stp_trace_single(ax1, x, iterations_remote, time_remote, x_label, y_label, 'r', 'Migratable Job (remote)', [time_pred[i]+time_own[i] for i in range(0,ranks-1)])

  acc_tmp = [time_pred[i]+time_own[i]+time_remote[i] for i in range(0,ranks-1)]
  acc = 0
  for n in range(0,ranks-1):
      acc = acc+acc_tmp[n]

  legend = ax1.legend(loc='upper left')
  ax1.set_ylim(0,2200000)
  #ax1.set_title("accumulated "+str(acc))

  plt.savefig(dirname+"stp_time_"+str(timestep)+".pdf")
  import tikzplotlib
  tikzplotlib.save(dirname+"stp_time_"+str(timestep)+".tex")

