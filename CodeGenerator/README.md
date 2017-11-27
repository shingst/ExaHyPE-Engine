Prerequisites
=============

Python
------

Python 3.3 or more is required.

Dependencies
------------

The CodeGenerator uses the template engine Jinja2 (http://jinja.pocoo.org/)

If it is not provided by your python installation you can instead use the source 
of jinja2 directly either by

1) Running the bash script: `./CodeGenerator/importJinjaAndUseSource.sh`

or by following these steps

1) Clone the source from the git repository to your desired install directory <my-path>: 
`git clone https://github.com/pallets/jinja.git <my-path>/jinja`

2) In ExaHyPE-Engine/CodeGenerator, create a symbolic link to <my-path>/jinja.
		
3) Modify TemplatingUtils.py to use the local version of Jinja2

From: 
```
isJinjaAvailableAsPackage=True
```

To 
```
isJinjaAvailableAsPackage=False
```

You may need to adapt the path of the sys.import
```
sys.path.insert(0, 'jinja')
```

Dependencies of jinja2
----------------------

jinja2 depends on the python module markupsafe. If it is not available on your 
system, get it manually via:
```
cd <my-path>
git clone https://github.com/pallets/markupsafe.git markupsafe-install
```
Finally, create a symbolic link to <my-path>/markupsafe-install/markupsafe
in ExaHyPE-Engine/Codegenerator via:
```
ln -s <my-path>/markupsafe-install/markupsafe
```


Paths
-----

Every path is relative to the root of the project (inside the ExaHyPE-Engine directory).

The CodeGenerator assumes the following:

* it is located in ``CodeGenerator/``
* the internal ExaHyPE is at ``ExaHyPE/``

If this is not the case, you may need to edit

* the configuration parameters of ``Toolkit/src/eu/exahype/CodeGeneratorHelper.java``
* the configuration parameters of ``CodeGenerator/Driver.py``

The generated code will be put accordingly to the ``pathToOptKernel`` argument starting from the internal ExaHyPe, by default in a subdirectory of ``ExaHyPE/kernels/aderdg/optimised/``.


Codegenerator
=============

To access the help: ``python3 CodeGenerator/Drivers.py -h``

Usage and arguments
-------------------

positional arguments:
*  pathToOptKernel    desired relative path to the generated code (../ExaHyPE/ as root)
*  solverName         name of the user-solver
*  numberOfVariables  the number of quantities
*  order              the order of the approximation polynomial
*  dimension          number of dimensions you want to simulate
*  numerics           linear or nonlinear
*  architecture       the microarchitecture of the target device
*  pathToLibxsmm      where to find your local copy of code generator back end 'https://github.com/hfp/libxsmm'

optional arguments:
*  -h, --help         show this help message and exit
*  --deepProfiling    enable deep-rpofiling (use only with profiler enable)
*  --useFlux          enable flux
*  --useNCP           enable non conservative product
*  --useSource        enable source terms
*  --noTimeAveraging  disable time averaging in the spacetimepredictor (less memory usage, more computation)


Example: ``python3 CodeGenerator/Driver.py kernels/aderdg/optimised/test Euler::MyEulerSolver 5 6 2 nonlinear hsw Libxsmm --useFlux``


Data format and padding
-----------------------

The Codegenerator may use padding when producing architecture specific code, it may also change the index order

Using the C index order convention with index in italic being absent in dim 2


| Array | tempArray Name | Generic | Optimised | Note |
| ----- | -------------- | ------- | --------- | ---- | 
| luh | | _nDof_, nDof, nDof, nData | _nDof_, nDof, nDof, nData | unchanged for compatibility purpose|
| lQhbnd | PEANO data | 2*nDim, _nDof_, nDof, **nData** | 2*nDim, _nDof_, nDof, **nDataPad** | 2*nDim face on the square/cube |
| lFhbnd | PEANO data | 2*nDim, _nDof_, nDof, **nVar** | 2*nDim, _nDof_, nDof, **nVarPad** | 2*nDim face on the square/cube |
| lQhi | tempUnknowns | _nDof_, nDof, nDof, **nData** | _nDof_, nDof, nDof, **nDataPad** | |
| lFhi | tempFluxUnknowns | nDim+1, _nDof_, nDof, nDof, **nVar** | nDim+1, _nDof_, nDof, nDof, **nVarPad** | lFhi has nDim+1 blocks == each directions + source |
| LQi **NL** | tempSpaceTimeUnknowns[0]  | nDof, _nDof_, nDof, nDof, **nData** | _nDof_, nDof, nDof, nDof, **nDataPad** | nonlinear case |
| LQi **Lin** | tempSpaceTimeUnknowns[0] | nDof+1, _nDof_, nDof, nDof, **nData** | nDof+1, _nDof_, nDof, nDof, **nDataPad** | linear case |
| LFi **NL** | tempSpaceTimeFluxUnknowns[0]  | nDof, _nDof_, nDof, nDof, nDim+1, **nVar** | nDim+1, nDof, _nDof_, nDof, nDof, **nVarPad** | nonlinear case, +1 for source. Move dimension coordinate to mimick lFhi blocks |
| LFi **Lin** | tempSpaceTimeFluxUnknowns[0] | 2*nDim+1, nDof, _nDof_, nDof, nDof, **nVar** | 2*nDim+1, nDof, _nDof_, nDof, nDof, **nVarPad** | linear case, +1 for source. |
| rhs | tempSpaceTimeUnknowns[1] | _nDof_, nDof, nDof, nDof, **nData** | _nDof_, nDof, nDof, nDof, **nDataPad** | nonlinear case |
| gradQ | tempSpaceTimeFluxUnknowns[1] | _nDof_, nDof, nDof, nDof, nDim, **nVar** | _nDof_, nDof, nDof, nDof, nDim, **nVarPad** | Not used if no NCP |
