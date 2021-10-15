#include <iostream>
#include "GridSolver.h"

using namespace std;

int main(int argc, char *argv[])
{
  if (argc == 2)
  {
    GridSolver gridSolver(argv[1]);
    gridSolver.printSolution();
  }
  else
  {
    cout << "Please provide a file name as the only command line argument for this program." << endl;
  }

#ifdef _WIN32
  system("PAUSE");
#endif
  
  return 0;
}