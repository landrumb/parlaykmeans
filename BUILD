licenses(["notice"])

cc_binary(
  name="txt_to_fbin",
  srcs=["txt_to_fbin.cpp"],
  linkopts=["-pthread"],
  deps = [
    "@parlaylib//parlay:parallel",
    "@parlaylib//parlay:sequence",
    "@parlaylib//parlay:primitives",
    "@parlaylib//parlay:io",

  ],
)

cc_binary(
    name = "test_main",
    srcs = ["test.cpp"],
    linkopts = ["-pthread"],
    deps = [
      "@parlaylib//parlay:parallel",
      "@parlaylib//parlay:sequence",
      "@parlaylib//parlay:primitives",
      "kmeans_headers" #remove when I finish testing NSG_Dist
    ],
)


#What's with this binary library separation?
cc_binary(
  name="kmeans_test_run",
  srcs=["kmeans.cpp"],
  linkopts=["-pthread"],
  #makes it known that include is an include library
  copts = ["-Iinclude"],
  deps = [
    ":kmeans_test",
  ],
 
)


cc_library(
  name="kmeans_test",
  srcs=["kmeans.cpp"],
  linkopts=["-pthread"],
  #makes it known that include is an include library
  copts = ["-Iinclude"],
  deps = [
    "kmeans_headers",

  ],
 
)

#for benching initializers
cc_binary(
  name="init_test_run",
  srcs=["bench_initializers.cpp"],
  linkopts=["-pthread"],
  #makes it known that include is an include library
  copts = ["-Iinclude"],
  deps = [
    ":kmeans_headers",

  ],
 
)


#for benching update step 
cc_binary(
  name="update_bench_run",
  srcs=["bench_update.cpp"],
  linkopts=["-pthread"],
  #makes it known that include is an include library
  copts = ["-Iinclude"],
  deps = [
    ":kmeans_headers",

  ],
 
)

#for benching assign step
cc_binary(
  name="assign_bench_run",
  srcs=["benchmarks/assign_bench.cpp"],
  linkopts=["-pthread"],
  #makes it known that include is an include library
  copts = ["-Iinclude"],
  deps = [
    ":kmeans_headers",

  ],
 
)
#for benching assign step
cc_binary(
  name="group_by_bench_run",
  srcs=["benchmarks/bench_group_by.cpp"],
  linkopts=["-pthread"],
  #makes it known that include is an include library
  copts = ["-Iinclude"],
  deps = [
    ":kmeans_headers",

  ],
 
)


cc_binary(
  name="debug",
  srcs=["debug.cpp"],
  linkopts=["-pthread"],
  #makes it known that include is an include library
  copts = ["-Iinclude"],
  deps = [
    ":kmeans_headers",

  ],
 
)


#comparing a yinyang and naive run
# TODO: this does not seem to build
cc_binary(
  name="kmeans_test_yy_naive",
  srcs=["kmeans_old.cpp"],
  linkopts=["-pthread"],
  #makes it known that include is an include library
  copts = ["-Iinclude"],
  deps = [
    ":kmeans_headers",
  ],
 
)



#comparing a yinyang and naive run
#but the yy version that doesn't local filter
# TODO: this does not seem to build
cc_binary(
  name="kmeans_test_yy_simp_naive",
  srcs=["kmeans_old.cpp"],
  linkopts=["-pthread"],
  #makes it known that include is an include library
  copts = ["-Iinclude"],
  deps = [
    ":kmeans_headers",
  ],
 
)


#Testing the tester
cc_test(
  name = "kmeans_google_test",
  size = "small",
  srcs = ["tests/extra_tests.cpp"],
 
  deps = ["@googletest//:gtest_main",
    "kmeans_headers"],
  linkopts=["-pthread"],
  #makes it known that include is an include library
  copts = ["-Iinclude"],
  
)

cc_test(
  name="naive_gtest",
  size="small",
  srcs=["test_n2.cpp"],
  deps=["@googletest//:gtest_main",
    "kmeans_headers"],
  linkopts=["-pthread"],
  copts = ["-Iinclude"],
)

cc_test(
  name="distance_gtest",
  size="small",
  srcs=["tests/distanceTests.cpp"],
  deps=["@googletest//:gtest_main",
    "kmeans_headers"],
  linkopts=["-pthread"],
  copts = ["-Iinclude"],
)

cc_library(
  name="kmeans_headers",
  hdrs=["include/utils/parse_files.h",
"include/lazy.h",
"include/utils/NSGDist.h",
"include/initialization.h",
"include/naive.h",
"include/utils/parse_command_line.h",
"include/utils/kmeans_bench.h",
"include/utils/threadlocal.h",
"include/yinyang_simp.h",
"include/utils/union_find.h",
"include/yy_structs.h",
"include/yy_compute_centers.h",
"include/quantized.h",
"include/nisk_kmeans.h",
"include/lsh.h",
"include/lsh_quantized.h",
"include/yy_improved/yy.h",
"include/yy_improved/yy_compute_centers.h",
"include/yy_improved/yy_structs.h",
"include/pq_kmeans.h",
"include/naive2.h"
],

linkopts=["-pthread"],
#makes it known that include is an include library
copts = ["-Iinclude"],
deps = [
  "@parlaylib//parlay:parallel",
  "@parlaylib//parlay:primitives",
  "@parlaylib//parlay:sequence",
  "@parlaylib//parlay:slice",
  "@parlaylib//parlay:io",
],

)

#Testing the tester
cc_test(
  name = "yy_google_test",
  size = "small",
  srcs = ["yy_tests.cpp"],
  deps = ["@googletest//:gtest_main",
    "kmeans_headers"],
  linkopts=["-pthread"],
  #makes it known that include is an include library
  copts = ["-Iinclude"],
  
)

#Testing the tester
cc_test(
  name = "nisk_google_test",
  size = "small",
  srcs = ["nisk_tests.cpp"],
  deps = ["@googletest//:gtest_main",
    "kmeans_headers"],
  linkopts=["-pthread"],
  #makes it known that include is an include library
  copts = ["-Iinclude"],
  
)



#Testing the tester
cc_test(
  name = "kmeans_plus_plus_google_test",
  size = "small",
  srcs = ["kmeans_plus_plus_tests.cpp"],
  deps = ["@googletest//:gtest_main",
    "kmeans_headers"],
  linkopts=["-pthread"],
  #makes it known that include is an include library
  copts = ["-Iinclude"],
  
)

#get a text version of a data file
cc_binary(
    name = "bin_to_csv",
    srcs = ["include/utils/convert_bin_to_csv.cpp"],
    linkopts = ["-pthread"],
    deps = [
     "parse_headers"
    ],
)



cc_library(
  name="parse_headers",
  hdrs = ["include/utils/parse_files.h",],

linkopts=["-pthread"],
#makes it known that include is an include library
copts = ["-Iinclude"],

  deps = [
  "@parlaylib//parlay:parallel",
  "@parlaylib//parlay:primitives",
  "@parlaylib//parlay:sequence",
  "@parlaylib//parlay:slice",
  "@parlaylib//parlay:io",
],
)

cc_binary(
  name="int_sort_practice_run",
  srcs=["int_sort_practice.cpp",],
  linkopts=["-pthread"],
  deps=["kmeans_headers"],

)
