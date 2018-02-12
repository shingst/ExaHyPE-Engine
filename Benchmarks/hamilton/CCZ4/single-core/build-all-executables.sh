export MODE='Release'
export COMPILER='Intel'
export SHAREDMEM='None'
export DISTRIBUTEDMEM='None'
export USE_IPO='off'

CORES=28

DRIVER='../../../CodeGenerator/Driver.py'

FOLDER=single-core
SPEC=${FOLDER}/CCZ4-no-output.exahype
APP=ExaHyPE-CCZ4
NEW_APP=ExaHyPE-CCZ4
#ARCH=knl
ARCH=hsw

# A1 -- A2.3
# pre: remove  -OX -ip -fma from the Makefile
declare -a CFLAGS=(' -O2 -no-vec -no-fma -no-ip -no-simd ' ' -O2 -fma -ip ' ' -O3 -fma -ip ' ' -O2 -fma -ip -DUSE_IPO ' ' -O3 -fma -ip -DUSE_IPO ' ' -O2 -fma -ip -DUSE_IPO ' ' -O3 -fma -ip -DUSE_IPO ' ' -O2 -fma -ip -DUSE_IPO ' ' -O3 -fma -ip -DUSE_IPO ') # last two rows for optimised

declare -a kernels=('generic' 'generic' 'generic' 'generic' 'generic' 'optimised' 'optimised' 'optimised' 'optimised')

declare -a suffixes=('O2-novec' 'O2-vec' 'O3-vec' 'O2-vec-ipo' 'O3-vec-ipo' 'opt-O2-vec-ipo-no-libxsmm' 'opt-O3-vec-ipo-no-libxsmm' 'opt-02-vec-ipo-libxsmm' 'opt-03-vec-ipo-libxsmm')

for arch in noarch $ARCH
do 
  sed -i -r 's,architecture(\s*)const(\s*)=(\s*).+,architecture\1const\2=\3'$arch',' $SPEC
  
  for i in {0..8}
  do   
      if (( i>=7  )); then
        sed -i -r 's,(\s*)'useLibxsmm'(\s*)\:(\s*)(True|False),\1'useLibxsmm'\2:\3True,' $DRIVER
      else
        sed -i -r 's,(\s*)'useLibxsmm'(\s*)\:(\s*)(True|False),\1'useLibxsmm'\2:\3False,' $DRIVER
      fi

      suffix=${suffixes[i]} 
      kernel=${kernels[i]}   

      echo ''
      echo $arch'-'${CFLAGS[i]}': '$suffix
      echo ''
 
    if [ "$arch" == "$ARCH" ] || [ "$kernel" == "generic" ]; then
      sed -i -r 's, optimisation(\s*)const(\s*)=(\s*)(generic|optimised), optimisation\1const\2=\3'$kernel',' $SPEC
      export COMPILER_CFLAGS=${CFLAGS[i]}
      make clean 
      
      for order in 3 6 9
      do
        rm *.o cipofiles.mk ffiles.mk cfiles.mk
        sed -i -r 's,order(\s*)const(\s*)=(\s*).+,order\1const\2=\3'$order',' $SPEC
        ${FOLDER}/configure-no-output.sh
        cat $SPEC

        make -j$cores
        mv $APP ${NEW_APP}-no-output-${arch}-$suffix-p${order}

      done
    fi   
  done
done
