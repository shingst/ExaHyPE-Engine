# ExaHyPE ApplicationExamples

This directory holds a large number of applications (implemented and tested PDEs).
The following list shall give an overview, with documentation:

```
ApplicationExamples
├── AcousticWave  (Seismology)
│   ├── acousticwaveCPP
│   ├── acousticwaveF90
│   └── PML_AcousticWaves_3D
├── DIM  (Diffusive Interface Method (Trento))
│   ├── DIM
│   ├── DIM_3DLADG
│   ├── DIM_FV
│   └── DIM_LimitingADERDG
├── Elastic   (Seismology)
│   ├── Elastic2D_Topo
│   ├── Elastic3D
│   ├── elasticwaveCPP
│   ├── elasticwaveF90
│   ├── elasticwaves2d
│   ├── ElasticWaves2D_nonlinearkernel
│   ├── elasticwaves3d
│   ├── ElasticWaves3D_nonlinearkernel
│   ├── ElasticWavesK2D
│   ├── ElasticWavesK2D_Probe
│   ├── ElasticWavesK3D_Probe
│   ├── ElasticWavesLOH1
│   ├── Elastodynamics2D
│   ├── Elastodynamics3D
│   └── LOH1
│       ├── LOH12D
│       ├── LOH1_2D
│       ├── LOH1_3D
│       └── LOH1_3D_probe
├── Euler              (The famous EulerFlow basic test)
│   ├── Euler          (Limiting application)
│   ├── Euler_ADERDG   (pure ADERDG application)
│   └── Euler_FV       (pure FV application)
├── EulerFlow          (Another collection of EulerFlow variants)
│   ├── doc
│   ├── EulerFlow
│   ├── EulerFlow_2ADERDG
│   ├── EulerFlow_f90
│   ├── EulerFlow_FV
│   ├── EulerFlow_H5
│   ├── EulerFlow_LimitingADERDG
│   ├── EulerFlow_NC
│   ├── EulerFlow_OptimisedKernel
│   └── Euler_FV_png
├── Experiments       (All kind of experimental and smallish dummy testing applications)
│   ├── Advection
│   ├── CellularAutomata
│   ├── DummyLimitingSolver
│   ├── DummySolver
│   ├── GridDemonstrator
│   ├── LinearAdvection
│   ├── PotentialHydro
│   ├── scalarwaveF90
│   └── SyntheticalTrouble
├── GRMHD           (The General-Relativistic Magnetohydrodynamics PDEs)
│   ├── doc         (LaTeX documentation)
│   ├── GRMHD       (Trento Fortran implementation)
│   ├── GRMHD_cpp   (C++ implementation)
│   ├── GRMHD_f90   (Trento Fortran with Fortran kernels, will not work any more)
│   ├── Specfiles
│   └── SVEC        (The C++ tensor code used for GRMHD_cpp)
├── Linear   (Seismology)
│   ├── ElasticWave2D
│   ├── ElasticWave2D_Topo
│   ├── ElasticWave3D
│   ├── ElasticWave3D_easi
│   ├── ElasticWave3D_Rupture
│   ├── ElasticWave3D_Topo
│   ├── LinearWave2D
│   ├── LinearWave2D_Topo
│   ├── LinearWave3D
│   ├── LinearWave3D_PML
│   └── LinearWave3D_Topo
├── Linear_Topo (Seismology)
├── SRHD           (The Special relativistic hydrodynamics, currently not maintained)
│   ├── SRHD
│   ├── SRHD_FV
│   └── SRHD_LimitingADERDG
├── SRMHD          (The Special relativistic magnetohydrodynamics, no more maintained)
│   ├── SRMHD
│   ├── SRMHD_cpp
│   ├── SRMHD_f90
│   ├── SRMHD_FV
│   ├── SRMHD_LimitingADERDG
│   └── SRMHD_TwoKernels
└── SWE            (The Shallow Water equations)
    ├── SWE_ADERDG
    └── SWE_FV
```