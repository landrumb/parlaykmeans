#!bin/bash

#file for benching initializers 

bazel build init_test_run

for multi in {10,10,200,400,500,1000}
do
echo "multi: c$multi"


./bazel-bin/init_test_run -k $multi -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c ForgyRandom -o "./experiments/experiment17/output/"

done


for multi in {2000,3125,4000,5000}
do
echo "multi: c$multi"

mya=$((1000 * $multi)) #no spaces!!

echo "a: $mya"

./bazel-bin/init_test_run -k $multi -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -c ForgyRandom -o "./experiments/experiment17/output/"

done