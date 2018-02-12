directory=multicore

APP=ExaHyPE-ElasticWaveEquation3D
NEW_APP=ExaHyPE-Elastic3D
SPEC=$directory/Elastic3D-no-output.exahype

# save original file
cp $SPEC ${SPEC}_tmp

for mode in None TBB OMP
do
  make clean
  export SHAREDMEM=$mode
  echo "SHAREDMEM=$SHAREDMEM"
  #read -p "press any key..."

  for p in 3 5 7 9
  do 
    rm *.o
    sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$p',' $SPEC
    cat $SPEC
    $directory/configure-no-output.sh
    make -j24 && \
    mv $APP $APP-p$p-$SHAREDMEM-$COMPILER
  done
done

# restore original file
mv ${SPEC}_tmp $SPEC
