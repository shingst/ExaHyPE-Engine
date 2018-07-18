Run:

1) Symlink ExaHyPE-Engine/Benchmarks/python/sweep.py into project folder (or create an alias)

2) Call ./sweep.py sweeps/dim_3dladgom.ini build

3) Symlink the geometry *.bin file into the sweeps/dim_3dladgom/build/ folder

4) Call ./sweep.py sweeps/dim_3dladgom.ini scripts

5) Call ./sweep.py sweeps/dim_3dladgom.ini submit

Parse output files:

1) Call ./sweep.py sweeps/dim_3dladgom.ini parseAdapters [--compress] # stuff in [...] is optional

2) Call ./sweep.py sweeps/dim_3dladgom.ini parseTimeStepTimes         # in order to get time per ADER-DG time step
