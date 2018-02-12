
export MODE="Release"
export COMPILER="Intel"
export SHAREDMEM="None"
export DISTRIBUTEDMEM="None"
# export USE_IPO="off"

CORES=28

FOLDER=a2.0
SPEC=${FOLDER}/CCZ4-no-output.exahype
APP=ExaHyPE-CCZ4
#ARCH=knl
ARCH=hsw

# A1
# pre: remove -O from the Makefile

for opt in "optimised" "generic"
do
	for ipo in "off" "on"
	do
	  export USE_IPO=$ipo   
	  for i in 0 1 2
	  do  
	    export COMPILER_CFLAGS="TODO"
	    make clean
	 
	    for order in 3 6 9
	    do
	      rm *.o cipofiles.mk ffiles.mk cfiles.mk Fortran/*.o
	      sed -i -r 's,order(\s*)const(\s*)=(\s*).+,order\1const\2=\3'$order',' $SPEC
	      sed -i -r 's,optimisation(\s*)const(\s*)=(\s*).+,optimisation\1const\2=\3'$opt',' $SPEC

	      
	      ${FOLDER}/configure-no-output.sh
	      cat $SPEC

	      # O2 no-vec no-fma
	      make -j$cores
	      mv $APP $APP-$opt-ipo-$ipo-p$order

	    done   
	  done
	done 
done

