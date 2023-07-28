#!bin/bash

#full experiment comparing naive and yy
echo "Make sure you move any previous experiments out of the output folder first."
echo "Run experiment? [Y/n]"
read input 
if [ $input = "Y" ] 
then 
  echo "Running the experiment6" 
  bazel build kmeans_test_run
  echo "Finishing compiling"

  # ./bazel-bin/kmeans_test_run -k 10 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D short -m 20 -two yes -csv_log true -csv_log_file_name experiment4/output/test_lazy_naive_10.csv -csv_log_file_name2 experiment4/output/test_lazy_yy_10.csv

 ./bazel-bin/kmeans_test_run -k 10000 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D short -m 20 -two yes -csv_log true -csv_log_file_name experiment7/output/test_lazy_quant_10000.csv 

  # ./bazel-bin/kmeans_test_run -k 1000 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D short -m 20 -two yes -csv_log true -csv_log_file_name experiment4/output/test_lazy_naive_1000.csv -csv_log_file_name2 experiment4/output/test_lazy_yy_1000.csv

  # ./bazel-bin/kmeans_test_run -k 10000 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D short -m 20 -two yes -csv_log true -csv_log_file_name experiment4/test_lazy_naive_10000.csv -csv_log_file_name2 experiment4/test_lazy_yy_10000.csv

  #  ./bazel-bin/kmeans_test_run -k 100000 -i ../my_data/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D short -m 20 -two yes -csv_log true -csv_log_file_name experiment4/test_lazy_naive_100000.csv -csv_log_file_name2 experiment4/test_lazy_yy_100000.csv

  # python3 include/utils/comparative_time_graphs.py < experiment3/comparative_plot_testing_input.txt

  # python3 include/utils/graph_bench2.py < experiment3/naive_iter_input.txt

  # python3 include/utils/graph_bench2.py < experiment3/yy_iter_input.txt

else
  echo "Not running the experiment6"

fi