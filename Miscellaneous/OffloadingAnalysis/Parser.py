import os.path

def timeSum( path ):
    mySum=0
    if(os.path.exists(path)):
        data = open(path,'r').read()
        values=data.split('\n')
        i=0
        while(i<len(values)):
            try:
                mySum+=int(values[i])
            except ValueError:
                print(values[i])
            i=i+1
    else:
        mySum=-1
    return str(mySum)



prefix="../../../AstroApplications/CCZ4/TraceOutput/"
rank=1
region=0
#regions=["exahype_solvers_ADERDGSolver_PredictionJob_run","iter_counter"]
#regions=["exahype_solvers_ADERDGSolver_MigratablePredictionJob_run","exahype_solvers_ADERDGSolver_PredictionJob_run","iter_counter"]
regions=["exahype_solvers_ADERDGSolver_MigratablePredictionJob_run","exahype_solvers_ADERDGSolver_PredictionJob_run","iter_counter","exahype_solvers_ADERDGSolver_PredictionJob_iterations","exahype_solvers_ADERDGSolver_MigratableJob_iterations"]

while(region<len(regions)):
    rank=1
    print()
    print(regions[region])
    print("Thread 0")
    while(rank<28):
        f= prefix+regions[region]+"_rank_"+str(rank)+"_"+"0.txt"        
        print(timeSum(f),end='')
        rank=rank+1
    rank=1
    print()
    print("Thread 1")
    while(rank<28):
        f= prefix+regions[region]+"_rank_"+str(rank)+"_"+"1.txt"
        print(timeSum(f),end='')
        rank=rank+1
    region=region+1
    print()
