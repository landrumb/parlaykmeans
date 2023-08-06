#!bin/bash

#full experiment comparing naive and yy
echo "Make sure you move any previous experiments out of the output folder first."
echo "Run experiment? [Y/n]"
read input 
if [ $input = "Y" ] 
then 
  echo "Running the experiment16" 
  bazel build kmeans_test_run
  echo "Finishing compiling"

  for myk in {1000,2000,3125,4000,5000} 
  do

  echo "multi: c$myk"


  ./bazel-bin/kmeans_test_run -k $myk -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D short -m 5 -two yes -csv_log true -csv_log_file_name "experiments/experiment16/output/" -csv_log_file_name2 "experiments/experiment16/output/"

  done


else
  echo "Not running the experiment16"

fi