/* Copyright (C) University of Utah - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Michael Matheny <michaelmathen@gmail.com>, 2015-2016
 */
#include <boost/format.hpp>
#include <sstream>
#include "Point.hpp"
namespace anomaly {
  Point::Point(){}

  Point::Point(double input_x, double input_y) : x(input_x),
						 y(input_y),
						 label(false)
  {}

  Point::Point(double input_x, double input_y, bool label) : x(input_x),
							     y(input_y),
							     label(label)
  {}

  void Point::set(double input_x, double input_y){
    x = input_x;
    y = input_y;
  }

  void Point::setAnomaly(){
    label = true;
  }

  string Point::str() {
    ostringstream os;
    os << *this;
    return os.str();
  }
  
  void Point::setNormal(){
    label = false;
  }

  void Point::setX(double input_x){
    x = input_x;
  }

  void Point::setY(double input_y){
    y = input_y;
  }

  bool Point::getLabel() const {
    return label;
  }

  bool Point::operator==(const Point& pOther) const {
    return sameLoc(pOther) &&
      (label == pOther.label);
  }

  ostream& operator<<(ostream& os, const Point& pt) {
    
    #ifdef DEBUG_POINT
      return os << (boost::format("Point(%.1f, %.1f, %s %d)") %
		    pt.x %
		    pt.y %
		    pt.text %
		    pt.label);
    #else
      return os << (boost::format("Point(%.1f, %.1f, %d)") %
		    pt.x %
		    pt.y %
		    pt.label);
    #endif
  }

  void Point::print() const {
    cout << "(" << x << ", " << y << ")" << endl;
  }

}
