#!bin/bash

#full experiment comparing naive and yy
echo "Make sure you move any previous experiments out of the output folder first."
echo "Run experiment? [Y/n]"
read input 
if [ $input = "Y" ] 
then 
  echo "Running the experiment5" 
  bazel build kmeans_test_run
  echo "Finishing compiling"

  ./bazel-bin/kmeans_test_run -k 10000 -i /ssd1/data/bigann/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D short -m 20 -two yes -csv_log true -csv_log_file_name experiment5/output/test_lazy_naive_10000.csv -csv_log_file_name2 experiment5/output/test_lazy_yy_10000.csv


else
  echo "Not running the experiment5"

fi