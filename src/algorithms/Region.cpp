/* Copyright (C) University of Utah - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Michael Matheny <michaelmathen@gmail.com>, 2015-2016
 */
#include <limits>

#include <boost/format.hpp>
#include <sstream>
#include "Region.hpp"

namespace anomaly {

  Disk::Disk(double a, double b, double r) :
    a(a),
    b(b),
    radius(r)
  {}

  Disk::Disk(Point const& pt) : radius(0) {
    pt.getLoc(a, b);
  }
  
  Disk::Disk(Point const& pt1, Point const& pt2) {
    radius = pt1.distanceWith(pt2);
    pt1.getLoc(a, b);
  }

  
  Disk::Disk(Point const& pt1, Point const& pt2, Point const& pt3) : Region(){
    if (pt1.sameLoc(pt2) && pt1.sameLoc(pt3)) {
      a = pt1.getX();
      b = pt1.getY();
      radius = 0;
    } else if (pt1.sameLoc(pt2) || pt2.sameLoc(pt3)) {
      solveCircle2(pt1, pt3, a, b, radius);
    } else if (pt1.sameLoc(pt3)) {
      solveCircle2(pt1, pt2, a, b, radius);
    }  else {
      solveCircle3(pt1, pt2, pt3, a, b, radius);
      //If they all lie on a line then a and b will be null
      //In this case the disk is a line
      if (std::isnan(a) || std::isnan(b)) {
	a = numeric_limits<double>::infinity();
	b = numeric_limits<double>::infinity();
	radius = numeric_limits<double>::infinity();
	x0 = pt1.getX();
	y0 = pt1.getY();
	u = pt1.getX() - pt2.getX();
	v = pt1.getY() - pt2.getY();
      }
    }
  }

  string Disk::str() {
    ostringstream os;
    os << *this;
    return os.str();
  }

  string Rectangle::str() {
    ostringstream os;
    os << *this;
    return os.str();
  }

  Rectangle::Rectangle(double x_min, double x_max, double y_min, double y_max)
    : x_max(x_max), x_min(x_min), y_min(y_min), y_max(y_max)
  {}

  bool Rectangle::inside(Point const& pt1) const{
    return pt1.getX() <= x_max &&
      x_min <=pt1.getX() &&
      pt1.getY() <= y_max &&
      y_min <=pt1.getY();
  }
  bool Disk::inside(Point const& pt1) const {
    return inDisk(a, b, radius, pt1);
  }
  
  void Region::setNumAnomalies(int in, int out) {
    numAnomalies = in;
    totalAnom = out;
  }
  void Region::setNumPoints(int in, int out) {
    numPoints = in;
    totalPoints = out;
  }

  int Region::getNumPoints() const {
    return numPoints;
  }

  int Region::getNumAnomalies() const {
    return numAnomalies;
  }

  int Region::getTotalPoints() const {
    return totalPoints;
  }

  int Region::getTotalAnomalies() const {
    return totalAnom;
  }
  ostream& operator<<(ostream& os, const Disk& dsk) {
    return os << (boost::format("Disk, %.1f, %.1f, %.1f, %d, %d %d %d") %
		  dsk.a %
		  dsk.b %
		  dsk.radius %
		  dsk.numAnomalies %
		  dsk.totalAnom %
		  dsk.numPoints %
		  dsk.totalPoints );
  }

  ostream& operator<<(ostream& os, const Rectangle& rect) {
    return os << (boost::format("Rectangle((%.1f, %.1f), (%.1f, %.1f), %d, %d, %d, %d)") %
		  rect.x_max %
		  rect.y_max %
		  rect.x_min %
		  rect.y_min %
		  rect.numAnomalies %
		  rect.totalAnom %
		  rect.numPoints %
		  rect.totalPoints);
  }

}
