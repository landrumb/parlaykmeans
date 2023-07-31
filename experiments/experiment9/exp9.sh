#!bin/bash

#full experiment comparing naive and yy
echo "Make sure you move any previous experiments out of the output folder first."
echo "Run experiment 9? [Y/n]"
read input 
if [ $input = "Y" ] 
then 
  echo "Running the experiment9" 
  bazel build kmeans_test_run
  echo "Finishing compiling"

  ./bazel-bin/kmeans_test_run -k 10 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D short -m 5 -two yes -csv_log true -csv_log_file_name experiment9/output/test_lazy_naive_fern_10.csv 

  PARLAY_NUM_THREADS=1 ./bazel-bin/kmeans_test_run -k 10 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D short -m 5 -two yes -csv_log true -csv_log_file_name experiment9/output/test_lazy_naive_fern_10_onethread.csv 


else
  echo "Not running the experiment9"

fi