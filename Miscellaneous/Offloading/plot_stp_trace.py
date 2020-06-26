#1/usr/bin/env python

import parse_stp_trace as pt
import sys
import matplotlib.pyplot as plt

chart_name="Euler Scenario Time Spent in STPs per Rank"
x_label="Ranks"
y_label="Elapsed time in seconds"

basename1 = "exahype_solvers_ADERDGSolver_PredictionJob__run_iterations_rank_"
basename2 = "exahype__solvers_ADERDGSolver_OwnMigratableJob__run_iterations_rank_"
basename3 = "exahype__solvers_ADERDGSolver_AlienMigratableJob__run_iterations_rank_"


def plot_stp_trace_single(ax,x, iterations, time, x_label, y_label, color, bottom):
  ax.set_xlabel(x_label)
  ax.set_ylabel(y_label)
  #ax.title(chart_name)
  ax.bar(x,time,0.8,label="Picard Iterations",bottom=bottom, color=color)
  #plt.bar(ranks,values_fts,0.8,Color="Red",bottom=values_stp,label="FTS")
  #plt.get_xaxis().set_ticks([0,2,4,6,8,10,12,14,16,18,20,22,24,26])
  #ax.set_ylim(0,600)

if __name__=="__main__":
  dirname = sys.argv[1]
  ranks = int(sys.argv[2])
  threads = int(sys.argv[3])  
 
  f, ax1 = plt.subplots(1)

  (x, iterations_pred, time_pred) = pt.parseDir(dirname, ranks, threads, basename1) 
  plot_stp_trace_single(ax1, x, iterations_pred, time_pred, x_label, y_label, 'b',0)
   
 # (x, iterations_own, time_own) = pt.parseDir(dirname, ranks, threads, basename2) 
 # plot_stp_trace_single(ax1, x, iterations_own, time_own, x_label, y_label, 'g',time_pred)

 # (x, iterations_remote, time_remote) = pt.parseDir(dirname, ranks, threads, basename3) 
 # plot_stp_trace_single(ax1, x, iterations_remote, time_remote, x_label, y_label, 'r',  [time_pred[i]+time_own[i] for i in range(0,ranks-1)])
  
  plt.show()
  
  legend = ax.legend(loc='upper left')

