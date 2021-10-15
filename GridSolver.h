#ifndef GRID_SOLVER_H
#define GRID_SOLVER_H

#include <list>
#include <string>
#include <unordered_set>
#include <vector>
#include "Point.h"

#define EPSILON 0.12

using namespace std;

class GridSolver
{
public:
  GridSolver(string fileName);
  ~GridSolver();

  void printSolution();

  enum Comparison
  {
    ComparisonLessThan,
    ComparisonEqual,
    ComparisonGreaterThan
  };

private:
  void           processGridFile(string fileName);
  void           buildBalancedKdTree(vector<Point*> &points, Point* &currentNode, int subtreeStartIdx, int subtreeEndIdx, int depth);
  void           quickSortPoints(vector<Point*> &points, int subtreeLowIdx, int subtreeHighIdx, int depth);
  int            quickSortPartition(vector<Point*> &points, int low, int high, int depth);
  void           swap(vector<Point*> &points, int idxA, int idxB);
  void           processGridSolution();
  vector<Point*> findNearestNeighbors(Point* &target);
  void           findNearestNeighbors(Point* &parent, Point* &target, int depth, vector<Point*> &nearestNeighbors);
  Comparison     compare(Point* &p0, Point* &p1, Point* &target);
  Point*         getCandidateBestNearestNeighbor(vector<Point*> &nearestNeighbors);
  double         calculateDistance(Point* &p0, Point* &p1);
  double         calculateAngle(Point* &p0, Point* &p1);
  void           deleteTreeNodes(Point* root);

  Point*             mTreeRoot;
  Point*             mTopLeftPoint;
  unordered_set<int> mVisitedPoints;
};

#endif