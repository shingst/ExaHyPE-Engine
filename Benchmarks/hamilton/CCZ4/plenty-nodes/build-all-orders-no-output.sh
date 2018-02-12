directory=plenty-nodes

exe=ExaHyPE-CCZ4
spec=$directory/CCZ4-no-output.exahype

# save original file
cp $spec ${spec}_tmp

for mode in None TBB OMP
do
  make clean
  export SHAREDMEM=$mode
  echo "SHAREDMEM=$SHAREDMEM"
  #read -p "press any key..."

  for p in 3 5 7
  do 
    rm *.o
    sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$p',' $spec
    cat $spec
    $directory/configure-no-output.sh
    make -j24 && \
    mv $exe $exe-p$p-$SHAREDMEM-$COMPILER
  done
done

# restore original file
mv ${spec}_tmp $spec
