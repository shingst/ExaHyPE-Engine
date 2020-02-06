def timeSum( path ):
    mySum=0
    data = open(path,'r').read()
    values=data.split('\n')
    i=0
    while(i<len(values)):
        try:
            mySum+=int(values[i])
        except ValueError:
            print(values[i])
        i=i+1
    return mySum



prefix="../../ExaHyPE-Max/ApplicationExamples/Euler/Euler_ADERDG/TraceOutput/"
rank=1
region=0
regions=["exahype_solvers_ADERDGSolver_PredictionJob_run","iter_counter"]
#regions=["exahype_solvers_ADERDGSolver_MigratablePredictionJob_run","exahype_solvers_ADERDGSolver_PredictionJob_run","iter_counter"]

while(region<len(regions)):
    rank=1
    print(regions[region])
    print("Thread 0")
    while(rank<28):
        f= prefix+regions[region]+"_rank_"+str(rank)+"_"+"0.txt"        
        print(timeSum(f),end='')
        rank=rank+1
    rank=1
    print("Thread 1")
    while(rank<28):
        f= prefix+regions[region]+"_rank_"+str(rank)+"_"+"1.txt"
        print(timeSum(f),end='')
        rank=rank+1
    region=region+1
