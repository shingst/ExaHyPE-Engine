export MODE="Release"
export COMPILER="Intel"
export SHAREDMEM="None"
export DISTRIBUTEDMEM="None"
export USE_IPO="off"

CORES=28

DRIVER="../../../CodeGenerator/Driver.py"

FOLDER=single-core
SPEC=${FOLDER}/Elastic3D-no-output.exahype
APP=ExaHyPE-ElasticWaveEquation3D
#ARCH=knl
ARCH=hsw

# A1 -- A2.3
# pre: remove  -OX -ip -fma from the Makefile
# pre: disable LIBXSMM in CodeGenerator/Driver.py
declare -a CFLAGS=(" -O2 -no-vec -no-fma -no-ip -no-simd" \ 
                   " -O2 -fma -ip " \
                   " -O3 -fma -ip " \ 
                   " -O3 -fma -ip -DUSE_IPO " \ 
                   " -O3 -fma -ip -no-simd"  \ 
                   " -O3 -fma -ip -DUSE_IPO -no-simd"  \ 
                   " -O3 -fma -ip -DUSE_IPO "  \ 
                   " -O3 -fma -ip -DUSE_IPO " )           # last four rows for optimised

declare -a KERNELS=( "generic" "generic" "generic" "generic" "optimised" "optimised" "optimised" "optimised" )

declare -a suffixes=("O2-novec" "O2-vec" "O3-vec" "O3-vec-ipo" "opt-O3-vec-no-ipo-no-simd-no-libxsmm" "opt-O3-vec-ipo-no-simd-no-libsxmm" "opt-O3-vec-ipo-simd-no-libxsmm" "opt-03-vec-ipo-simd-libxsmm" )

for arch in noarch $ARCH
do 
  sed -i -r 's,architecture(\s*)const(\s*)=(\s*).+,architecture\1const\2=\3'$arch',' $SPEC
  
  for i in {0..3}
  do  
    if (( i==7 )); then
      sed -i -r 's,(\s*)"useLibxsmm"(\s*)\:(\s*)(True|False),\1"useLibxsmm"\2:\3True,' $DRIVER
    else
      sed -i -r 's,(\s*)"useLibxsmm"(\s*)\:(\s*)(True|False),\1"useLibxsmm"\2:\3False,' $DRIVER
    fi

    kernel=${KERNELS[i]}   
    suffix=${suffixes[i]} 

    sed -i -r 's, optimisation(\s*)const(\s*)=(\s*)(generic|optimised), optimisation\1const\2=\3'$kernel',' $SPEC
    export COMPILER_CFLAGS="${CFLAGS[i]}"
    make clean
 
    for order in 3 6 9
    do
      rm *.o cipofiles.mk ffiles.mk cfiles.mk
      sed -i -r 's,order(\s*)const(\s*)=(\s*).+,order\1const\2=\3'$order',' $SPEC
      
      ${FOLDER}/configure-no-output.sh
      cat $SPEC

      # O2 no-vec no-fma
      make -j$cores
      mv $APP $APP-no-output-${arch}-$suffix-p${order}

    done   
  done
done
