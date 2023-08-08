
bazel build update_bench_run



for multi in {10,100,1000,10000,100000}
do
echo "multi: c$multi"

./bazel-bin/update_bench_run -k $multi -i /ssd1/anndata/spacev1b_base.i8bin -f bin -t uint8 -D Euclidean -c Standard -v add_overn -o "./experiments/experiment18/output/"

done


