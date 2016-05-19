#ifndef ANOMALY_REGION
#define ANOMALY_REGION
#include <iostream>

#include "Utility.hpp"
#include "Point.hpp"
namespace anomaly {


  class Disk;
  
  class Region {
  protected:
    int numAnomalies = 0;
    int numPoints = 0;
    int totalAnom = 0;
    int totalPoints = 0;
    bool usePrecalculated = false;
    double stat = 0;
  public:
    Region(){}
    virtual void setNumAnomalies(int in, int out);
    virtual void setNumPoints(int in, int out);
    virtual int getNumAnomalies() const;
    virtual int getNumPoints() const;
    virtual bool inside(Point const& pt1) const = 0;
    virtual int getTotalAnomalies() const;
    virtual int getTotalPoints() const;
    virtual void useKStat(double d) {
      usePrecalculated = true;
      stat = d;
    }
    virtual double statistic() {
      if (usePrecalculated) {
	return stat;
      } else {
	return kulldorff(getNumPoints(),
			 getTotalPoints(),
			 getNumAnomalies(),
			 getTotalAnomalies());
      }
    }
    friend ostream& operator<<(ostream& os, const Disk& pt);
  };
  
  class Disk : public Region {
    double a;
    double b;
    double radius;
    double u = 0, v = 0, x0 = 0, y0 = 0;
  public:
    Disk(){}
    Disk(double a, double b, double r);
    Disk(Point const& pt);
    Disk(Point const& pt1, Point const& pt2);
    Disk(Point const& pt1, Point const& pt2, Point const& pt3);
    bool inside(Point const& pt1) const;
    friend ostream& operator<<(ostream& os, const Disk& pt);
    double getX() { return a; }
    double getY() { return b; }
    double getRadius() { return radius; }

    void setX(double x) { a = x; }
    void setY(double y) { b = y; }
    void setRadius(double r) { radius = r; }
    
    string str();
  };

  class Rectangle : public Region {
    double x_min, y_min;
    double x_max, y_max;
  public:
    double getLX() { return x_min; }
    double getLY() { return y_min; }
    double getUX() { return x_max; }
    double getUY() { return y_max; }

    void setLX(double x) { x_min = x; }
    void setLY(double y) { y_min = y; }
    void setUX(double x) { x_max = x; }
    void setUY(double y) { y_max = y; }

    Rectangle(){}
    Rectangle(double x_min, double x_max, double y_min, double y_max);
    string str();
    bool inside(Point const& pt1) const;
    friend ostream& operator<<(ostream& os, const Rectangle& pt);
  };

  ostream& operator<<(ostream& os, const Rectangle& pt);

  ostream& operator<<(ostream& os, const Disk& pt);

}
#endif
