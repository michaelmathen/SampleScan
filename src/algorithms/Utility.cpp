#include "Utility.hpp"

namespace anomaly {

  void solveCircle3(Point const& pt1,
		    Point const& pt2,
		    Point const& pt3,
		    double &a,
		    double &b,
		    double &r) {
    double x1, x2, x3, y1, y2, y3;
    pt1.getLoc(x1, y1);
    pt2.getLoc(x2, y2);
    pt3.getLoc(x3, y3);
    // Setup a matrix equation of the form Ax = b

    // A
    double
      a11 = x2 - x1, a12 = y2 - y1,
      a21 = x2 - x3, a22 = y2 - y3;
    // b
    double
      b1 = (y2 * y2 + x2 * x2 - y1 * y1 - x1 * x1) / 2,
      b2 = (y2 * y2 + x2 * x2 - y3 * y3 - x3 * x3) / 2;

    double detA = a11 * a22 - a12 * a21;
    // inverse of A
    double
      ai11 =  a22 / detA, ai12 = -a12 / detA,
      ai21 = -a21 / detA, ai22 =  a11 / detA;

    // A^{-1} b = x
    a = ai11 * b1 + ai12 * b2;
    b = ai21 * b1 + ai22 * b2;
    r = sqrt((x1 - a) * (x1 - a) + (y1 - b) * (y1 - b));
  }


  bool colinear(Point const& pt1,
		Point const& pt2,
		Point const& pt3){
    double x1, x2, x3, y1, y2, y3;
    pt1.getLoc(x1, y1);
    pt2.getLoc(x2, y2);
    pt3.getLoc(x3, y3);

    double
      a11 = x2 - x1, a12 = y2 - y1,
      a21 = x2 - x3, a22 = y2 - y3;
    return (a11 * a22 - a12 * a21 == 0);
  }

  bool onLineSegment(Point const& pt1,
		     Point const& pt2,
		     Point const& pt3) {
    if (colinear(pt1, pt2, pt3)) {
      if (pt1.sameLoc(pt2))
	return pt1.sameLoc(pt3);
      //Now we know that the point is on the same line
      if (pt1.getX() != pt2.getX()) {
	double theta = (pt1.getY() - pt3.getX()) / (pt1.getX() - pt2.getX());
	return (theta <= 1) && (theta >= 0);
      } else {
	double theta = (pt1.getY() - pt3.getY()) / (pt1.getY() - pt2.getY());
	return (theta <= 1) && (theta >= 0);
      }
    } else {
      return false;
    }
  }


  void solveCircle3(Point const& pt1,
		    Point const& pt2,
		    Point const& pt3,
		    double &a,
		    double &b) {
    double x1, x2, x3, y1, y2, y3;
    pt1.getLoc(x1, y1);
    pt2.getLoc(x2, y2);
    pt3.getLoc(x3, y3);
    // Setup a matrix equation of the form Ax = b

    // A
    double
      a11 = x2 - x1, a12 = y2 - y1,
      a21 = x2 - x3, a22 = y2 - y3;
    // b
    double
      b1 = (y2 * y2 + x2 * x2 - y1 * y1 - x1 * x1) / 2,
      b2 = (y2 * y2 + x2 * x2 - y3 * y3 - x3 * x3) / 2;

    double detA = a11 * a22 - a12 * a21;
    // inverse of A
    double
      ai11 =  a22 / detA, ai12 = -a12 / detA,
      ai21 = -a21 / detA, ai22 =  a11 / detA;

    // A^{-1} b = x
    a = ai11 * b1 + ai12 * b2;
    b = ai21 * b1 + ai22 * b2;
  }

  

  void solveCircle2(Point const& pt1,
		    Point const& pt2,
		    double &a,
		    double &b,
		    double &r) {
    double x1, x2, y1, y2;
    pt1.getLoc(x1, y1);
    pt2.getLoc(x2, y2);
  
    r = pt1.distanceWith(pt2) / 2.0;
    a = (x1 + x2) / 2;
    b = (y1 + y2) / 2;
  }


  bool inDisk(double a, double b, double r, Point const& pt) {
    double x, y;
    pt.getLoc(x, y);
    return (x - a) * (x - a) + (y - b) * (y - b) <= r * r;
  }



  void findPerpVect(Point const& p1, Point const& p2, double* u, double* v) {
    double x1, x2, y1, y2;
    p1.getLoc(x1, y1);
    p2.getLoc(x2, y2);
    *u = y2 - y1;
    *v = x1 - x2;
  }


  void findOrthoVector(double at1, double at2,
		       double &x1, double &y1) {
  
    // Generate a random point uniformly around a disk centered at the midpoint.
    // between the two edge points.
    double r = std::sqrt(at1 * at1 + at2 * at2) / 2;
    std::default_random_engine generator;
    std::normal_distribution<double> dist(0, r);
    
    // We generate a random point inside of a disk of radius r. Think of this as a
    // random vector emanating from the midpoint between our two reference points.
    // we then check to make sure that this vector isn't dependent on the vector
    // between the reference points (close to 0 probability of this happenning)
    // We can then orthogonalize against our vector between the reference points.
    double proj, norm, x, y;
    do {

      // We want to ensure that this newly generated vector does not
      // lie on the vector between these two points because this would imply it is dependent.
      x = dist(generator);
      y = dist(generator);

      // Project the vector <x1, y1> onto the vector lieing between the two reference points.
      proj = std::abs((x * at1  + y * at2) / std::sqrt(at1 * at1 + at2 * at2));
      norm = std::sqrt(x * x + y * y);
    } while (std::abs((proj - norm) / norm) < 1e-6);

    // Now orthoganilize against the midpoint 
    x = x - x * at1 / std::sqrt(at1 * at1 + at2 * at2);
    y = y - y * at2 / std::sqrt(at1 * at1 + at2 * at2);
      
    // and normalize
    double tmpNorm = std::sqrt(x * x + y * y);
    x1 = x / tmpNorm;
    y1 = y / tmpNorm;

  }



  void parsePtFile(string const& fname, vector<Point>& pts){
    ifstream ptFile;
    ptFile.open (fname);
    string lineBuff;
    do {
      vector<string> token_parts;
      getline(ptFile, lineBuff);
      boost::split(token_parts, lineBuff, boost::is_any_of(","));
      if (token_parts.size() == 2) {
	pts.push_back(Point(stod(token_parts[0]),
			    stod(token_parts[1])));
      } else if (token_parts.size() == 3) {
	pts.push_back(Point(stod(token_parts[0]),
			    stod(token_parts[1]),
			    stoi(token_parts[2])));
      } else if (token_parts.size() == 0 || token_parts[0] == "") {
	break;
      } else {
	throw "Something weird happened";
      }
    } while (!ptFile.eof());
    ptFile.close();
}
}

