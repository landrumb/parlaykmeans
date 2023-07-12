#include <stdio.h>
#include <fstream>
#include <ostream>
#include <iostream>
#include <cmath> 
#include "parlay/parallel.h"

int main() {
  std::cout << "hello" << std::endl;
  std::cout << parlay::num_workers() << std::endl;

  // sequence<size_t> test;
  // test.push_back(5);
  // test.push_back(22);
  // auto dist = parlay::map(test,[&](size_t j))
  std::cout << "min of " << std::min(5,7) << std::endl;
  size_t a = 0;
  size_t b = 1;
  std::cout << "min of " << std::min(a,b) << std::endl;
  return 1;
}
