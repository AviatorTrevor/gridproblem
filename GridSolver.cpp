#include <fstream>
#include <iostream>
#include <math.h>
#include "GridSolver.h"
#include "Point.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
GridSolver::GridSolver(string fileName) :
  mTreeRoot(nullptr),
  mTopLeftPoint(nullptr)
{
  processGridFile(fileName);
  processGridSolution(mTreeRoot);
}

///////////////////////////////////////////////////////////////////////////////
GridSolver::~GridSolver()
{
  deleteTreeNodes(mTreeRoot);
}

///////////////////////////////////////////////////////////////////////////////
// recursive function, iterates top-to-bottom in the tree and deletes nodes to
// "properly" clean up the memory we allocated on the heap.
///////////////////////////////////////////////////////////////////////////////
void GridSolver::deleteTreeNodes(Point* root)
{
  if (!root->isNull())
  {
    deleteTreeNodes(*root->getLeftChild());
    deleteTreeNodes(*root->getRightChild());
    delete root;
  }
}

///////////////////////////////////////////////////////////////////////////////
// This function reads the file, creates the "Point" objects, and inserts
// them into a k-d tree
///////////////////////////////////////////////////////////////////////////////
void GridSolver::processGridFile(string fileName)
{
  ifstream myFile;
  string myLine;
  int idCounter = 0;
  vector<Point*> points;
  vector<Point*>::iterator it;
  myFile.open(fileName, ios::in);
  if (myFile.is_open())
  {
    while (!myFile.eof())
    {
      getline(myFile, myLine);
      if (myLine.length() > 2) //has to be a length greater than 2, because you have 2 numbers and a comma separator. This helps when the last line of the file is empty
      {
        for (unsigned int i = 0; i < myLine.length(); i++)
        {
          if (myLine[i] == ',')
          {
            //Extract the coordinates; convert them from strings to doubles
            double x = stod(myLine.substr(0, i)); //NOTE: I am not validating the input. For the purpose of this exercise, I assume the input is a valid number
            double y = stod(myLine.substr(i + 1, myLine.length() - i + 1));

            //Create a new Point
            Point* p = new Point(x, y, idCounter++);

            //insertion-sort; sorting by x-coordinate to prepare it for making a balanced kd-tree
            for (it = points.begin(); it != points.end(); it++)
            {
              if (p->getX() < (*it)->getX())
              {
                points.insert(it, p);
                break;
              }
            }
            if (it == points.end())
              points.push_back(p);
            break;
          }
        }
      }
    }
  }
  else
  {
    cout << "Error opening file called \"" << fileName << "\"" << endl;
  }

  buildBalancedKdTree(points, mTreeRoot, 0, points.size() - 1, 0);
}

///////////////////////////////////////////////////////////////////////////////
// This function assumes the list of points has already been sorted by x-coordinate
// For larger grids, a balanced tree will help with run-time efficiency
///////////////////////////////////////////////////////////////////////////////
void GridSolver::buildBalancedKdTree(vector<Point*> &points, Point* &currentNode, int subtreeStartIdx, int subtreeEndIdx, int depth)
{
  int median = subtreeStartIdx + (subtreeEndIdx - subtreeStartIdx + 1) / 2;
  currentNode = points[median]; //set the current node to the median

  if (subtreeStartIdx == subtreeEndIdx)
    return;
  
  //Sort the left-branch and continue with recursively creating the balanced-tree
  if (median != subtreeStartIdx)
  {
    quickSortPoints(points, subtreeStartIdx, median - 1, depth + 1);
    buildBalancedKdTree(points, *(currentNode->getLeftChild()), subtreeStartIdx, median - 1, depth + 1);
  }

  //Sort the right-branch and continue with recursively creating the balanced-tree
  if (median != subtreeEndIdx)
  {
    quickSortPoints(points, median + 1, subtreeEndIdx, depth + 1);
    buildBalancedKdTree(points, *(currentNode->getRightChild()), median + 1, subtreeEndIdx, depth + 1);
  }
}

///////////////////////////////////////////////////////////////////////////////
// Implements a quick-sort program for a sub-section of the vector of points.
// It sorts those points by the dimension corresponding to the depth in the tree
///////////////////////////////////////////////////////////////////////////////
void GridSolver::quickSortPoints(vector<Point*> &points, int subtreeLowIdx, int subtreeHighIdx, int depth)
{
  if (subtreeLowIdx < subtreeHighIdx)
  {
    int pi = quickSortPartition(points, subtreeLowIdx, subtreeHighIdx, depth);

    quickSortPoints(points, subtreeLowIdx, pi - 1, depth);
    quickSortPoints(points, pi + 1, subtreeHighIdx, depth);
  }
}

///////////////////////////////////////////////////////////////////////////////
// This is called by the quickSortPoints function to aid in the
// quick-sort algorithm
///////////////////////////////////////////////////////////////////////////////
int GridSolver::quickSortPartition(vector<Point*> &points, int low, int high, int depth)
{
  double pivot = points[high]->getCoordinate(depth % POINT_DIMENSIONS); //get either the x or the y coordinate based on tree depth
  int i = low - 1;

  for (int j = low; j < high; j++)
  {
    if (points[j]->getCoordinate(depth % POINT_DIMENSIONS) <= pivot)
    {
      i++;
      swap(points, i, j);
    }
  }
  swap(points, i + 1, high);
  return i + 1;
}

///////////////////////////////////////////////////////////////////////////////
// This function swaps two Point* pointers in a vector given 2 indexes
///////////////////////////////////////////////////////////////////////////////
void GridSolver::swap(vector<Point*> &points, int idxA, int idxB)
{
  Point* temp = points[idxA];
  points[idxA] = points[idxB];
  points[idxB] = temp;
}

///////////////////////////////////////////////////////////////////////////////
// After the data has been entered into the k-d tree, the constructor calls
// this function to solve the entire grid. "Solving" entails setting the
// neighbor pointers for each point and then noting which point is at the
// top-left corner of the grid
///////////////////////////////////////////////////////////////////////////////
void GridSolver::processGridSolution(Point* node)
{
  mVisitedPoints.insert(node->getId());
  vector<Point*> nearestNeighbors = findNearestNeighbors(node);

  cout << "Visited [" << node->getX() << "," << node->getY() << "]. Neighbors:  "; //TODO: debug, remove
  for (Point* neighbor : nearestNeighbors)
  {
    cout << "[" << neighbor->getX() << "," << neighbor->getY() << "], "; //TODO: debug, remove
  }
  cout << endl; //TODO: debug, remove
  for (Point* neighbor : nearestNeighbors)
  {
    if (mVisitedPoints.count(neighbor->getId()) == 0)
    {
      processGridSolution(neighbor);
    }
  }

  //TODO calculateAlpha() at the end. Take the mTopLeftNode and it's right-side neighbor
}

///////////////////////////////////////////////////////////////////////////////
// It's probably always easier to call this function and not the other
// "findNearestNeighbors" function, that way you get the right starting
// parameters for the parent and depth
///////////////////////////////////////////////////////////////////////////////
vector<Point*> GridSolver::findNearestNeighbors(Point* &target)
{
  vector<Point*> nearestNeighbors;
  findNearestNeighbors(mTreeRoot, target, 0, nearestNeighbors);
  return nearestNeighbors;
}

///////////////////////////////////////////////////////////////////////////////
// This function should only be called by the *other* "findNearestNeighbors"
// function or by itself. It just keeps it cleaner since the initial target
// and depth parameters should always be the same.
//
// If the depth is even-numbered, we compare the X-coordinate to determine
// if we should go to the left or right child next. Otherwise, we use
// the Y-coordinate to determine if we should go to the left or right child
// next. When we get to the bottom of the tree we will have found our
// candidate nearest neighbor point, but it's possible that once we get to
// the bottom node that our sibling or its children might actually be closer,
// so we need to check that as well if the vertical or horizontal distance
// to the parent node is less than the node we currently found as the best
// match.
///////////////////////////////////////////////////////////////////////////////
void GridSolver::findNearestNeighbors(Point* &parent, Point* &target, int depth, vector<Point*> &nearestNeighbors)
{
  if (parent->isNull())
    return;

  Point* nextBranch;
  Point* siblingBranch;
  int d = depth % POINT_DIMENSIONS;

  if (target->getCoordinate(d) < parent->getCoordinate(d))
  {
    nextBranch = *(parent->getLeftChild());
    siblingBranch = *(parent->getRightChild());
  }
  else
  {
    nextBranch = *(parent->getRightChild());
    siblingBranch = *(parent->getLeftChild());
  }

  int d = depth % POINT_DIMENSIONS;
  findNearestNeighbors(nextBranch, target, depth + 1, nearestNeighbors);
  if (nearestNeighbors.empty())
  {
    nearestNeighbors.push_back(parent);
  }
  else
  {
    Point* best = closest(nearestNeighbors.back(), parent, target);
    if (best->isNull())
    {
      nearestNeighbors.pop_back();
    }
  }

  //There's a chance we should traverse the sibling and its children. This is how we check if we should go that way
  double distanceToBestSquared = calculateDistanceSquared(nearestNeighbors.front(), target);
  double distancePrime = abs(target->getCoordinate(d) - parent->getCoordinate(d));

  /* If the perpendicular line to the sibling branch cutting plane (distancePrime) is <= to our curren best
   * distance, we should check that sibling branch as well. Because of issues with calculating distances
   * between 2 points using floating points, there is a chance we have 2 comparisons that are "equal"
   * as judged by a human, but not equal as judged by the computer because it is slightly off based on
   * how the number is stored. That's why we add "EPSILON" to the distanceToBestSquared variable (to
   * allow for a little wiggle room on what we consider "equal") */
  if (distancePrime * distancePrime <= distanceToBestSquared + EPSILON)
  {
    findNearestNeighbors(siblingBranch, target, depth + 1, nearestNeighbors);
    best = closest(temp, best, target);
  }
}

///////////////////////////////////////////////////////////////////////////////
// Determines which point is closest. If it's comparing against itself, we 
// return the "other" point so as to not include ourselves in the nearest-neighbor
// search results. If the two points being compared are of equal distance to
// the target point, then we return a null pointer to signify it was a tie
///////////////////////////////////////////////////////////////////////////////
Point* GridSolver::closest(Point* &p0, Point* &p1, Point* &target)
{
  if (p0->isNull()) return p1;
  if (p1->isNull()) return p0;

  double d0 = calculateDistanceSquared(p0, target);
  double d1 = calculateDistanceSquared(p1, target);

  /* if we are comparing the target to itself (distance is zero),
   * then we don't count it as the "best" or "closest" */
  if (abs(d0) < EPSILON)
    return p1;
  else if (abs(d1) < EPSILON)
    return p0;

  if (abs(d0 - d1) < EPSILON)
    return nullptr;

  if (d0 < d1) //TODO look at this function
    return p0;
  else
    return p1;
}

///////////////////////////////////////////////////////////////////////////////
// Calculates distance squared using Pythagorean theorem: a^2 + b^2 = c^2.
// In this case, result^2 = x^2 + y^2, where x is the x-axis distance
// between p0 and p1, and where y is the y-axis distance between p0 and p1.
// Instead of taking the squart root to calculate the actual distance, we'll
// leave the result in the form of distance-squared so as to save some CPU
// cycles. But the answer will still give a magnitude that gives us a way to
// judge in a relative fashion if one distance is greater than another.
//
// This function should always return a positive value.
///////////////////////////////////////////////////////////////////////////////
double GridSolver::calculateDistanceSquared(Point* &p0, Point* &p1)
{
  double xDistance = p0->getX() - p1->getX();
  double yDistance = p0->getY() - p1->getY();
  return pow(xDistance, 2) + pow(yDistance, 2);
}

///////////////////////////////////////////////////////////////////////////////
// Since the desired solution is to simply print the values in the order of:
// 1) All rows in order
// 2) All columns in order
// 3) the angle of the grid
// ...
// ...
// ... I will print the rows as I iterate through the grid, and store the
// column data for printing later.
///////////////////////////////////////////////////////////////////////////////
void GridSolver::printSolution()
{
  //TODO
}