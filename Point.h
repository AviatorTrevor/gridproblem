#ifndef POINT_H
#define POINT_H

#define POINT_DIMENSIONS 2 //currently only support 2D points

class Point
{
public:
  Point(double x, double y, int id);

  double  getX();
  double  getY();
  double  getCoordinate(int index);
  int     getId();

  void    setLeftNeighbor(Point* point);
  void    setTopNeighbor(Point* point);
  void    setRightNeighbor(Point* point);
  void    setBottomNeighbor(Point* point);

  bool    hasLeftNeighbor();
  bool    hasTopNeighbor();
  bool    hasRightNeighbor();
  bool    hasBottomNeighbor();

  void    addLeftChild(Point* point);
  void    addRightChild(Point* point);

  Point** getParent();
  Point** getLeftChild();
  Point** getRightChild();

  bool    isNull();
  bool    isEqual(Point* p);

private:
  void    setParent(Point* point);

  double mCoordinates[POINT_DIMENSIONS];
  int    mPointId;

  //Neighbors for the Grid
  Point* mLeftNeighbor;
  Point* mTopNeighbor;
  Point* mRightNeighbor;
  Point* mBottomNeighbor;

  //Parent and Children for the tree
  Point* mParent;
  Point* mLeftChild;
  Point* mRightChild;
};

#endif