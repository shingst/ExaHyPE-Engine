Workflow
========

Always work in the project directory (here: ApplicationExamples/DIM/DIM_3DLADGOM).

1) Build debug executable
-------------------------

Load modules and set debugging environment variables
```
source ddt/env-supermuc-debug.sh
```
Prepare the build
```
make clean # only necessary if environment or dimension was changed
ddt/configure.sh
```
Run make
```
make -j16
```

2) Debug the executable with ddt
--------------------------------
Start DDT
```
module load ddt
ddt &
```
Use DDT 

Follow this [guide](https://www.lrz.de/services/software/programmierung/ddt/).

* Application: 
  ```
  <exahype-root>/ApplicationExamples/DIM/DIM_3DLADGOM/ExaHyPE-DIM
  ```
* Arguments: 
  ```
  ddt/DIM-debug.exahype
  ```
* Working Directory: 
  ```
  <exahype-root>/ApplicationExamples/DIM/DIM_3DLADGOM 
  ```

As submission template file select:
```
 <exahype-root>/ApplicationExamples/DIM/DIM_3DLADGOM/ddt/loadleveler.qtf 
```
