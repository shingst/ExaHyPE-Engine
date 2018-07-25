Run:

1) Call ./sweep.py sweeps/dim_3dladgom.ini build

2) Symlink the geometry *.bin file into the sweeps/dim_3dladgom/build/ folder

3) Call ./sweep.py sweeps/dim_3dladgom.ini scripts

4) Call ./sweep.py sweeps/dim_3dladgom.ini submit

Parse output files:

1) Call ./sweep.py sweeps/dim_3dladgom.ini parseAdapters [--compress] # stuff in [...] is optional

2) Call ./sweep.py sweeps/dim_3dladgom.ini parseTimeStepTimes         # in order to get time per ADER-DG time step
