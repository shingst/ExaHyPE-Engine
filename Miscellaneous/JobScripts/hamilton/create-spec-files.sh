if [ "${1##*.}" != "exahype" ]
then
 echo "usage: ./create-spec-files.sh <file.exahype>"
 exit
fi

h=("0.038" "0.013" "0.00412" "0.00138") # >1/27, >1/81, >1/243, >1/729
T=("0.12" "0.04" "0.013" "0.0044") # ~0.12 ~0.12/3 ~0.12/9 ~0.12/27

for NRANKS_PER_NODE in 4 8 16
do

for i in 0 1 2 3
do 

sed -e 's,ranks-per-node:[0-9]\+,ranks-per-node:'"${NRANKS_PER_NODE}"',' \
    -e 's,end-time                 =.\+,end-time                 = '"${T[$i]}"',' \
    -e 's,maximum-mesh-size =.\+,maximum-mesh-size = '"${h[$i]}"',' \
    $1 > ${1%.exahype}-${NRANKS_PER_NODE}-$i.exahype # ${1%.exahype} removes the suffix ".exahype"

done
done
