#!bin/bash
bazel build kmeans_test_run

# ./bazel-bin/kmeans_test_run -k 50 -i Data/diff_europe.bin -f bin -t float -D short -m 10 -two yes 


# ./bazel-bin/kmeans_test_run -k 10 -i Data/europediff_csv.csv -f bin -t float -D short -m 10 -two yes 

#./bazel-bin/kmeans_test_run -k 20 -i ./Data/base.1B.u8bin.crop_nb_1000 -f bin -t uint8 -D short -m 5 -two yes 

# ./bazel-bin/kmeans_test_run -k 10 -i ./Data/base.1B.u8bin.crop_nb_1000 -f bin -t uint8 -D short -m 10 -two yes -csv_log true -csv_log_file_name plots/debugging/test_lazy_naive_10.csv

# ./bazel-bin/kmeans_test_run -k 10 -i ./Data/base.1B.u8bin.crop_nb_1000 -f bin -t uint8 -D short -m 10 -two yes -csv_log true -csv_log_file_name plots/debugging/test_lazy_naive_10.csv

#./bazel-bin/kmeans_test_run -k 20 -i /ssd1/anndata/bigann/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D short -m 5 -two yes 

./bazel-bin/kmeans_test_run -k 200 -i /ssd1/anndata/bigann/base.1B.u8bin.crop_nb_10000000 -f bin -t uint8 -D short -m 5 -two yes 
