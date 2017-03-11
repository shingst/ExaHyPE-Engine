import sys
import re
import os
from enum import Enum


Symbol = [ 
  "s", "o", ">", 
  "<", "^", "v" 
]
Colour = [
  "#ff0000", "#00ff00", "#0000ff",
  "#ffff00", "#ff00ff", "#00ffff"
]


class Analysis(Enum):
  Nop                = 0,
  NoSerialRuntimeYet = 1,
  SeemsNotToScale    = 2,
  MightScale         = 3,
  DoesNotScale       = 4,
  Scales             = 5


AdapterNames = []

def readAdapterNumbers( filename ):
  if not filename.endswith( "specification" ):
    print "WARNING: Peano specification files usually end with \"specification\"" 

  inputFile = open(filename, "r" )

  runIntoAdapter = False
  for line in inputFile:  
    if re.search( "adapter:", line ):
      runIntoAdapter = True
    elif re.search( "name:", line ) and runIntoAdapter:
      adapterName = line.split( "name:") [-1]
      AdapterNames.append( adapterName )
    else:
      runIntoAdapter = False
      

  print "finished to read specification file " + filename
  
  

def processMeasurement(adapter):
  searchPattern = "adapter=" + str(adapter) + ","

  searchPattern = "adapter-number=" + str(adapter) 
  line          = inputFile.readline()
  runtime       = 0
  while not line=="" and not re.search( searchPattern, line ):
    if re.search( "total-runtime=", line ):
      runtime = line.split("=")[1]
    line = inputFile.readline()
  
  htmlOverview.write( "<table border=\"1\">" );
  htmlOverview.write( 
    "<tr><td><b>Method</b></td>\
    <td><b>Maximum problem size</b></td>\
    \<td><b>Grain size</b></td>\
    <td><b>Search delta</b></td>\
    <td><b>Time serial</b></td>\
    <td><b>Time parallel</b>\
    <td><b>Accuracy serial</b></td>\
    <td><b>Accuracy parallel</b></td>\
    <td><b>Max. concurrency level</b></td>\
    <td><b>Predicted speedup</b></td>\
    <td><b>Speedup obtained before</b></td>\
    <td><b>Comments</b></td>\
    </tr>" );

  line = inputFile.readline()
  while not re.search( "end OracleForOnePhaseWithShrinkingGrainSize", line):
    methodTrace        = line.split( "=")[0]
    biggestProblemSize = line.split( "=")[1].split( "," )[0]
    grainSize          = line.split( "=")[1].split( "," )[1]
    searchDelta        = line.split( "=")[1].split( "," )[2]
    serialAccuracy     = line.split( "=")[1].split( "," )[3]
    parallellAccuracy  = line.split( "=")[1].split( "," )[4]
    speedupBefore      = line.split( "=")[1].split( "," )[5]
    timeSerial         = line.split( ",(" )[1].split( "," )[0]
    timeParallel       = line.split( ",(" )[2].split( "," )[0]
    
   
    analysis = Analysis.Nop
    if int(grainSize)==int(biggestProblemSize) and int(searchDelta)>0:
      analysis = Analysis.NoSerialRuntimeYet
    elif int(biggestProblemSize)==int(grainSize) and int(searchDelta)>0:
      analysis = Analysis.SeemsNotToScale
    elif int(biggestProblemSize)>int(grainSize) and int(searchDelta)>0:
      analysis = Analysis.MightScale
    elif int(biggestProblemSize)==int(grainSize) and int(searchDelta)==0:
      analysis = Analysis.DoesNotScale
    elif int(biggestProblemSize)>int(grainSize) and int(searchDelta)==0:
      analysis = Analysis.Scales


    if analysis==Analysis.NoSerialRuntimeYet and float(serialAccuracy)==0:
      colour = "White"
    elif analysis==Analysis.NoSerialRuntimeYet:
      colour = "Yellow"
    elif analysis==Analysis.SeemsNotToScale:
      colour = "Fuchsia"
    elif analysis==Analysis.MightScale:
      colour = "LightSkyBlue"
    elif analysis==Analysis.DoesNotScale:
      colour = "Red"
    elif analysis==Analysis.Scales:
      colour = "LightGreen"
    else:
      colour = "White"

    maxConcurrencyLevel = float(biggestProblemSize)/float(grainSize)
    if maxConcurrencyLevel<=1.0 and colour!="White":
      colour = "Grey"

    htmlOverview.write( "<tr>" );
    htmlOverview.write( "<td>" + methodTrace + "</td>" );
    htmlOverview.write( "<td>" + biggestProblemSize + "</td>" );
    htmlOverview.write( "<td bgcolor=\"" + colour + "\">" + grainSize + "</td>" );
    htmlOverview.write( "<td>" + searchDelta + "</td>" );
    htmlOverview.write( "<td>" + timeSerial + "</td>" );
    htmlOverview.write( "<td>" + timeParallel + "</td>" );
    htmlOverview.write( "<td>" + serialAccuracy + "</td>" );
    htmlOverview.write( "<td>" + parallellAccuracy + "</td>" );
    htmlOverview.write( "<td bgcolor=\"" + colour + "\">" + str(maxConcurrencyLevel) + "</td>" );

    if serialAccuracy!="nan" and float(serialAccuracy)>0.0:
      predicatedSpeedup = float(timeSerial)/float(timeParallel)
      if predicatedSpeedup<1.0:
        colour = "#ee8787"
      elif predicatedSpeedup>1.0:
        colour = "#87ee87"
      else:
        colour = "#878787"
      htmlOverview.write( "<td bgcolor=\"" + colour + "\">" + str( predicatedSpeedup ) + "</td>" );
    else:
      htmlOverview.write( "<td>no prediction</td>" );

    htmlOverview.write( "<td>" + speedupBefore + "</td>" );

    htmlOverview.write( "<td>");
    if int(searchDelta)!=0:
      htmlOverview.write( "Still searching. " );
    htmlOverview.write( str(analysis) );
    htmlOverview.write( ". " );
    if analysis==Analysis.SeemsNotToScale:
      htmlOverview.write( "Code might have found scaling setup but wants to re-validate serial runtime." );
    htmlOverview.write( "</td>" );    
    
    htmlOverview.write( "</tr>" );
    line = inputFile.readline()  

  htmlOverview.write( "</table>" );
      
  return;

#
#
#   main
# ========
#
#    
if (len(sys.argv)<3):
  print "usage: python ../postprocess-shrinking-grain-size-output.py Peano-specification-file output-files\n"
  print "\n"
  print "       You can hand over multiple output files in one rush\n"
  quit()
  
readAdapterNumbers( sys.argv[1] )

for i in range(2,len(sys.argv)):
  print "process data file " + sys.argv[i]
  htmlOverview = open( sys.argv[i] + ".html",  "w" )
  htmlOverview.write( "<h1>" + sys.argv[i] + "</h1>" );

  inputFile = open(sys.argv[i], "r" )

  line                  = inputFile.readline()
  totalNumberOfOracles  = int( line.split("=")[1] )
  line                  = inputFile.readline()
  oraclesForSteering    = int( line.split("=")[1] )

  htmlOverview.write( "Number of oracles=" + str(totalNumberOfOracles) + " (incl. oracles required for repository/algorithm steering)" );
  htmlOverview.write( "<br />" );
  htmlOverview.write( "Oracles required for repository steering=" + str(oraclesForSteering) );

  htmlOverview.write( "<h3>Table of content:</h3>" );
  htmlOverview.write( "<ul>" );
  for adapter in range(oraclesForSteering,totalNumberOfOracles):
    htmlOverview.write( "<li><a href=\"#adapter-" + str(adapter) + "\">Adapter " + str(adapter) + "</a></li>" );
  htmlOverview.write( "</ul>" );

  htmlOverview.write( "<p>Empty adapter sections imply that this adapter is not used by the code. The adapter order requals the adapter order in the specification file. Only those program phases that are actually used are also displayed.</p>" );
  htmlOverview.write( "<p>All individual timings are normalised by the number of entries handled, i.e. they do specify time per grid entity. We thus may not compare them directly to the total runtime.</p>" );

  for adapter in range(oraclesForSteering,totalNumberOfOracles):
    if adapter-oraclesForSteering>=len(AdapterNames):
      print "ERROR: No name available for adapter no " + str(adapter-oraclesForSteering)
      print "ERROR: Available adapters are " + str(AdapterNames)
      htmlOverview.write( "<h3 id=\"adapter-" + str(adapter) + "\">Adapter " + str(adapter) + "</h3>" );
    else:
      htmlOverview.write( "<h3 id=\"adapter-" + str(adapter) + "\">Adapter " + str(adapter) + ": " + AdapterNames[adapter-oraclesForSteering] + "</h3>" );
    processMeasurement(adapter)
    htmlOverview.write( "<p>The serial runtime is not reliable: It tracks only timinigs if the code fragment is not parallelised and the search for a grain size still is switched on. The longer your code runs, the more `invalid' this figure becomes. However, large entries indicate that you have heavy-weight serial part in your code.</p>" );
  