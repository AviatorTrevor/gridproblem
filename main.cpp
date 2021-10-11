#include <iostream>
#include "GridSolver.h"

using namespace std;

int main()
{
  //TODO: have the program argument be the text file input
  GridSolver gridSolver("grid_input.txt");
  gridSolver.printSolution();

#ifdef _WIN32
  system("PAUSE");
#endif
  

  //NOTE: have a hash map for "already-visited" points and do this recursively to expand out until all points have mapped their neighbors
  //NOTE: alpha angle can only be from 89.99999 to -89.999999
  return 0;
}