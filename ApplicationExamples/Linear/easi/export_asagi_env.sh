export COMPILER_CFLAGS="${PROJECT_CFLAGS} -DUSE_ASAGI"
export COMPILER_LFLAGS=""
export COMPILER_LFLAGS="${COMPILER_LFLAGS} -lasagi -lyaml-cpp -limpalajit -lnetcdf -lhdf5_hl -lhdf5 -lm -lz -lmpi -lcurl  -lnuma"
