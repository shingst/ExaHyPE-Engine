directory=plenty-nodes

APP=ExaHyPE-ElasticWaveEquation3D
NEW_APP=ExaHyPE-Elastic3D
SPEC=$directory/Elastic3D-output.exahype

# save originial file
cp $SPEC ${SPEC}_tmp

for m in 1 2
do
  if (( m == 1 )); then
    make clean
    export SHAREDMEM=TBB
  else
    make clean
    export SHAREDMEM=None
  fi

  echo "SHAREDMEM=$SHAREDMEM"
  #read -p "press any key..."

  for p in 9 7 5 3
  do 
    rm *.o
    sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$p',' $SPEC
    cat $SPEC
    $directory/configure-output.sh
    make -j24 && \
    mv $APP $NEW_APP-p$p-$SHAREDMEM-$COMPILER
  done
done

# restore original file
mv ${SPEC}_tmp $SPEC
