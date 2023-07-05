#include <stdio.h>
#include <fstream>
#include <ostream>
#include <iostream>

#include "parlay/parallel.h"

int main() {
  std::cout << "hello" << std::endl;
  std::cout << parlay::num_workers() << std::endl;
  return 1;
}
