/* Copyright (C) University of Utah - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Michael Matheny <michaelmathen@gmail.com>, 2015-2016
 */
#ifndef FRAMEWORK_POINT
#define FRAMEWORK_POINT
#include <string>
#include <cmath>
#include <cstdint>
#include <iostream>

namespace anomaly {
  using namespace std;
  class Point{
    double x, y;
    bool label;
    #ifdef DEBUG_POINT
    string text;
    #endif 
  public:
    Point();
    Point(double input_x, double input_y);
    Point(double input_x, double input_y, bool label);
    
    inline double distanceSquared(const Point& p1) const {
      double tmpx = x - p1.x;
      double tmpy = y - p1.y;
      return tmpx * tmpx + tmpy * tmpy;
    }
    
    inline double distanceWith(const Point& p1) const {
      double tmpx = x - p1.x;
      double tmpy = y - p1.y;
      return sqrt(tmpx * tmpx + tmpy * tmpy);
    }

    void set(double input_x, double input_y);
    void setAnomaly();
    void setNormal();
    void setLabel(bool label) {
      this->label = label;
    }
    void setX(double input_x);
    #ifdef DEBUG_POINT
    void setText(string text) {
      this->text = text;
    }
    #endif
    string str();
    void setY(double input_y);
    bool getLabel() const;
    bool sameLoc(const Point& pOther) const {
      return pOther.x == this->x && pOther.y == this->y;
    }
    
    inline void getLoc(double& x, double &y) const {
      x = this->x;
      y = this->y;

    }
    bool operator==(const Point& pOther) const;
    
    bool operator!=(const Point& pOther) const {
      return !(*this == pOther);
    }
    
    void print() const ;
    inline double getX() const {
      return x;
    }
    inline double getY() const {
      return y;
    }
    friend ostream& operator<<(ostream& os, const Point& pt);
  };

  ostream& operator<<(ostream& os, const Point& pt);
}
#endif

