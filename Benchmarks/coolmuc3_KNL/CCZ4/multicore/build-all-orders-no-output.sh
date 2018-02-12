directory=multicore

exe=ExaHyPE-CCZ4
spec=$directory/CCZ4-no-output.exahype

cp $spec ${spec}_tmp

for m in OMP TBB NONE
do
  make clean
  export SHAREDMEM=$m

  echo "SHAREDMEM=$SHAREDMEM"
  #read -p "press any key..."

  for p in 3 6 9
  do 
    rm *.o Fortran/*.o
    sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$p',' $spec
    cat $spec
    $directory/configure-no-output.sh
    make -j56 && \
    mv $exe $exe-p$p-$SHAREDMEM-$COMPILER
  done
done

mv ${spec}_tmp $spec
