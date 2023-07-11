#include <stdio.h>
#include <fstream>
#include <ostream>
#include <iostream>

#include "parlay/parallel.h"

int main() {
  std::cout << "hello" << std::endl;
  std::cout << parlay::num_workers() << std::endl;

  // sequence<size_t> test;
  // test.push_back(5);
  // test.push_back(22);
  // auto dist = parlay::map(test,[&](size_t j))
  return 1;
}
