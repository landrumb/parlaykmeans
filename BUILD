licenses(["notice"])

cc_binary(
    name = "test_main",
    srcs = ["test.cpp"],
    linkopts = ["-pthread"],
    deps = [
      "@parlaylib//parlay:parallel",
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
  hdrs = ["include/utils/parse_files.h",
  "include/lazy.h",
  "include/utils/NSGDist.h",
  "include/initialization.h",
  "include/naive.h",
  "include/utils/accumulator.h",
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


#What's with this binary library separation?
cc_binary(
  name="yinyang_debug_test",
  srcs=["include/yinyang_debugging.cpp"],
  linkopts=["-pthread"],
  #makes it known that include is an include library
  copts = ["-Iinclude"],
  deps = [
    ":yinyang_faithful",
  ],
 
)


#making sure that yinyang compiles to help with debugging
cc_library(
  name="yinyang_faithful",
  srcs=["include/yinyang_faithful.h"],
  hdrs=["include/NSGDist.h"],
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


