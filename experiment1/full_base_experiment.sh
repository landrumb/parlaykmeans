#!bin/bash

#full experiment comparing naive and yy

bazel build kmeans_test_run

./bazel-bin/kmeans_test_run -k 10 -i ./Data/base.1B.u8bin.crop_nb_1000 -f bin -t uint8 -D short -m 20 -two yes -csv_log true -csv_log_file_name experiment1/test_lazy_naive_10.csv -csv_log_file_name2 experiment1/test_lazy_yy_10.csv

./bazel-bin/kmeans_test_run -k 50 -i ./Data/base.1B.u8bin.crop_nb_1000 -f bin -t uint8 -D short -m 20 -two yes -csv_log true -csv_log_file_name experiment1/test_lazy_naive_50.csv -csv_log_file_name2 experiment1/test_lazy_yy_50.csv

./bazel-bin/kmeans_test_run -k 100 -i ./Data/base.1B.u8bin.crop_nb_1000 -f bin -t uint8 -D short -m 20 -two yes -csv_log true -csv_log_file_name experiment1/test_lazy_naive_100.csv -csv_log_file_name2 experiment1/test_lazy_yy_100.csv

python include/utils/comparative_time_graphs.py < experiment1/comparative_plot_testing_input.txt

python include/utils/graph_bench2.py < experiment1/naive_iternum_input.txt

python include/utils/graph_bench2.py < experiment1/yy_iternum_input.txt
