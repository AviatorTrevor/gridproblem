#ifndef GRID_SOLVER_H
#define GRID_SOLVER_H

#include <list>
#include <string>
#include <unordered_set>
#include <vector>
#include "Point.h"

#define EPSILON 0.01

using namespace std;

class GridSolver
{
public:
  GridSolver(string fileName);
  ~GridSolver();

  void printSolution();

private:
  void          processGridFile(string fileName);
  void          buildBalancedKdTree(vector<Point*> &points, Point* &currentNode, int subtreeStartIdx, int subtreeEndIdx, int depth);
  void          quickSortPoints(vector<Point*> &points, int subtreeLowIdx, int subtreeHighIdx, int depth);
  int           quickSortPartition(vector<Point*> &points, int low, int high, int depth);
  void          swap(vector<Point*> &points, int idxA, int idxB);
  void          insertPointIntoTree(Point* &parent, Point* newPoint, int depth);
  void          processGridSolution(Point* node);
  list<Point*>  findNearestNeighbors(Point* &target);
  void          findNearestNeighbors(Point* &parent, Point* &target, int depth, list<Point*> &nearestNeighbors);
  Point*        closest(Point* &p0, Point* &p1, Point* &target);
  double        calculateDistanceSquared(Point* &p0, Point* &p1);
  void          deleteTreeNodes(Point* root);

  Point*             mTreeRoot;
  Point*             mTopLeftPoint;
  unordered_set<int> mVisitedPoints;
};

#endif