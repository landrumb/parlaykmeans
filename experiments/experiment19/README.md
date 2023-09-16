In this experiment, experiment19, we will look at the runtime of eakmeans (Newling and Fleuret's library) across various 
values of k. We will examine the effects of running with 1 thread, no BLAS, 1 thread, yes BLAS, many (19) threads, no BLAS, and many (19) threads yes BLAS to try to understand to what extent these features speed up their code. Further, we will run differing versions of syin (syin-sn and syin-ns). We'll also collect the mse information and # of total distance calculations, to get time per distance calculation info.

We will do the same checks (with k=100 because I want these runs to finish) for sta/naive k-means/Lloyd's to help determine if this is an implementation problem (we did yy badly) or an us problem (we did both slowly).

Note that for yy sometimes I include the grouping simple k-means run, other times I don't (for brevity).

To continue this experiment, we run yy and naive with the vectorized distance class (which I estimate to be x4 faster), with all the threads, as well as without the vectorized distance class with 19 threads. 





Scratch:

bin/withblaskmeans -nc 1000 -din ../my_data/u8bin_1000000_txt -nth 19 --maxiter 5 -alg syin-ns
bin/blaslesskmeans -nc 5 -din examples/dense_200_5_header.txt -nth 19 --maxiter 5 -alg syin-ns
bin/blaslesskmeans -nc 1000 -din ../../my_data/u8bin_1000000_txt -nth 19 --maxiter 5 -alg syin-ns

** Use sn ! ** (for yy slower but more comparable to ours)

bin/blaslesskmeans -nc 100 -din ../../my_data/u8bin_1000000_txt -nth 1 --maxiter 5 -alg sta

bin/withblaskmeans -nc 100 -din ../my_data/u8bin_1000000_txt -nth 1 --maxiter 5 -alg sta

bazel build kmeans_test_run
./bazel-bin/kmeans_test_run -k 1000 -i /ssd1/anndata/bigann/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D Euclidean -m 5 -two yes 



bazel build kmeans_test_run
./bazel-bin/kmeans_test_run -k 1000 -i /ssd1/anndata/bigann/base.1B.u8bin.crop_nb_1000000 -f bin -t uint8 -D short -m 5 -two yes 
