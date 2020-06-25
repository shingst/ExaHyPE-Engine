#1/usr/bin/env python

import parse_stp_trace as pt
import sys
import matplotlib.pyplot as plt

chart_name="Euler Scenario Time Spent in STPs per Rank"
x_label="ranks"
y_label="time in Seconds"

def plot_stp_trace(x, iterations, time, chart_name, x_label, y_label):
  plt.xlabel(x_label)
  plt.ylabel(y_label)
  plt.title(chart_name)
  plt.bar(x,time,0.8,label="Picard Iterations")
  #plt.bar(ranks,values_fts,0.8,Color="Red",bottom=values_stp,label="FTS")
  legend = plt.legend(loc='upper left')
  #plt.get_xaxis().set_ticks([0,2,4,6,8,10,12,14,16,18,20,22,24,26])
  #ax.set_ylim(0,600)
  plt.show()

if __name__=="__main__":
  dirname = sys.argv[1]
  ranks = int(sys.argv[2])
  threads = int(sys.argv[3])  
 

  (x, iterations, time) = pt.parseDir(dirname, ranks, threads) 
  plot_stp_trace(x, iterations, time, chart_name, x_label, y_label)

