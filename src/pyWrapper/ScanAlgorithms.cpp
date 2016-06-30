/* Copyright (C) University of Utah - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Michael Matheny <michaelmathen@gmail.com>, 2015-2016
 */
#include <vector>

#include "algorithms/DiskScan.hpp"
#include "algorithms/Region.hpp"
#include "algorithms/RectScan.hpp"
#include "algorithms/Jeff_Algorithm.hpp"
#include "algorithms/Utility.hpp"

#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/stl_iterator.hpp>

boost::python::tuple perpVect(anomaly::Point& p1, anomaly::Point& p2) {
  double u, v;
  anomaly::findPerpVect(p1, p2, &u, &v);
  return boost::python::make_tuple(u, v);
}


boost::python::tuple disk3(anomaly::Point& p1, anomaly::Point& p2, anomaly::Point& p3) {
  double a, b, r;
  anomaly::solveCircle3(p1, p2, p3, a, b, r);
  return boost::python::make_tuple(a, b, r);
}

boost::python::tuple disk2(anomaly::Point& p1, anomaly::Point& p2, anomaly::Point& p3) {
  double a, b;
  anomaly::solveCircle3(p1, p2, p3, a, b);
  return boost::python::make_tuple(a, b);
}

template <typename T, typename ...Ss>
boost::python::object runAlg(boost::python::object const& ob, Ss ...args) {
  boost::python::stl_input_iterator<anomaly::Point> begin(ob), end;
  if (begin == end) {
    return boost::python::object();
  }
  T alg;
  return boost::python::object(alg(begin, end, args...));
}

boost::python::object runSatScan(boost::python::object const& ob) {
  boost::python::stl_input_iterator<anomaly::Point> begin(ob), end;
  if (begin == end) {
    return boost::python::object();
  }
  std::vector<anomaly::Point> pts(begin, end);
  anomaly::SatScan alg;
  return boost::python::object(alg(pts.begin(), pts.end()));
}


template <typename T>
boost::python::object runAlgSample(boost::python::object const& n,
				   boost::python::object const& s_a,
				   boost::python::object const& s_b,
				   double rho) {
  boost::python::stl_input_iterator<anomaly::Point> begin_n(n), end_n;
  boost::python::stl_input_iterator<anomaly::Point> begin_sa(s_a), end_sa;
  boost::python::stl_input_iterator<anomaly::Point> begin_sb(s_b), end_sb;
  if (begin_n == end_n ||
      begin_sa == end_sa ||
      begin_sb == end_sb) {
    return boost::python::object();
  }
  std::vector<anomaly::Point> n_pts(begin_n, end_n);
  std::vector<anomaly::Point> sa_pts(begin_sa, end_sa);
  std::vector<anomaly::Point> sb_pts(begin_sb, end_sb);  
  T alg;
  return boost::python::object(alg.run(n_pts.begin(), n_pts.end(),
				       sa_pts.begin(), sa_pts.end(),
				       sb_pts.begin(), sb_pts.end(),
				       rho));
}


boost::python::object runAlgSampleC(boost::python::object const& s_a,
				    boost::python::object const& s_b,
				    double rho) {
  boost::python::stl_input_iterator<anomaly::Point> begin_sa(s_a), end_sa;
  boost::python::stl_input_iterator<anomaly::Point> begin_sb(s_b), end_sb;
  if (begin_sa == end_sa ||
      begin_sb == end_sb) {
    return boost::python::object();
  }
  std::vector<anomaly::Point> sa_pts(begin_sa, end_sa);
  std::vector<anomaly::Point> sb_pts(begin_sb, end_sb);  
  anomaly::CenteredDisks alg;
  return boost::python::object(alg.run(sa_pts.begin(), sa_pts.end(), sb_pts.begin(), sb_pts.end(), rho));
}


boost::python::object runJeff(boost::python::object const& ob, double eps) {
  boost::python::stl_input_iterator<anomaly::Point> begin(ob), end;
  std::vector<anomaly::Point> newV(begin, end);
  if (begin == end) {
    return boost::python::object();
  }
  anomaly::JeffAlgo alg;
  return boost::python::object(alg(newV.begin(), newV.end(), eps));
}

inline anomaly::Disk make_disk3(anomaly::Point const& p1,
				anomaly::Point const& p2,
				anomaly::Point const& p3) {
  return anomaly::Disk(p1, p2, p3);
}

inline anomaly::Disk make_disk2(anomaly::Point const& p1,
				anomaly::Point const& p2) {
  return anomaly::Disk(p1, p2);
}

inline anomaly::Disk make_disk1(anomaly::Point const& p1) {
  return anomaly::Disk(p1);
}

BOOST_PYTHON_MODULE(eps_scan)
{
  using namespace boost::python;

  class_<anomaly::Region, boost::noncopyable>("Region", no_init)
    .def("numAnomalies", &anomaly::Region::getNumAnomalies)
    .def("numPoints", &anomaly::Region::getNumPoints)
    .def("setNumAnomalies", &anomaly::Region::setNumAnomalies)
    .def("setNumPoints", &anomaly::Region::setNumPoints)
    .def("totalPoints", &anomaly::Region::getTotalPoints)
    .def("totalAnomalies", &anomaly::Region::getTotalAnomalies)
    .def("statistic", &anomaly::Region::statistic)
    .def("inside", &anomaly::Region::inside)
    ;
  
  class_<anomaly::Disk, bases<anomaly::Region> >("Disk", init<double, double, double>())
    .add_property("a", &anomaly::Disk::getX, &anomaly::Disk::setX)
    .add_property("b", &anomaly::Disk::getY, &anomaly::Disk::setY)
    .add_property("radius", &anomaly::Disk::getRadius, &anomaly::Disk::setRadius)
    .def("__str__",  &anomaly::Disk::str)
    .def("__repr__",  &anomaly::Disk::str)
    ;

  class_<anomaly::Rectangle, bases<anomaly::Region> >("Rectangle", init<double, double, double, double>())
    .add_property("min_x", &anomaly::Rectangle::getLX, &anomaly::Rectangle::setLX)
    .add_property("min_y", &anomaly::Rectangle::getLY, &anomaly::Rectangle::setLY)
    .add_property("max_x", &anomaly::Rectangle::getUX, &anomaly::Rectangle::setUX)
    .add_property("max_y", &anomaly::Rectangle::getUY, &anomaly::Rectangle::setUY)
    .def("__str__",  &anomaly::Rectangle::str)
    .def("__repr__",  &anomaly::Rectangle::str)
    ;

  
  class_<anomaly::Point>("Point", init<double, double, bool>())
    .add_property("x", &anomaly::Point::getX, &anomaly::Point::setX)
    .add_property("y", &anomaly::Point::getY, &anomaly::Point::setY)
    .add_property("anomaly", &anomaly::Point::getLabel, &anomaly::Point::setLabel)
    .def("__str__",  &anomaly::Point::str)
    .def("__repr__",  &anomaly::Point::str)
    ;

  def("perpVect", &perpVect);
  def("onLineSegment", &anomaly::onLineSegment);
  def("colinear", &anomaly::colinear);
  def("kulldorff", &anomaly::kulldorff);
  
  def("disk2", &disk2);
  def("disk3", &disk3);

  def("makeDisk1", &make_disk1);
  def("makeDisk2", &make_disk2);
  def("makeDisk3", &make_disk3);
  
  def("netDisks", runAlg<anomaly::AllDisks, int, int, double>);
  //def("netRects", runAlg<anomaly::RectangleScan, int, int, double>);
  def("netRectsGrid", runAlg<anomaly::RectangleScanG, int, int, double>);
  def("netRectsEGrid", runAlg<anomaly::RectangleScanGH, int, int, double>);    
  def("netDisks2", runAlg<anomaly::MoreDisks, int, int, double>);
  def("netCDisks", runAlg<anomaly::CenterNetDisks, int, int, double>);
  def("cDisks", runAlg<anomaly::CenteredDisks, int, double>);
  def("scanRects", runJeff);
  def("satScan", runSatScan);
  
  def("netDisksSample", runAlgSample<anomaly::AllDisks>);
  //def("netRectsSample", runAlgSample<anomaly::RectangleScan>);
  def("netRectsGSample", runAlgSample<anomaly::RectangleScanG>);
  def("netRectsEGSample", runAlgSample<anomaly::RectangleScanGH>);
  def("netDisks2Sample", runAlgSample<anomaly::MoreDisks>);
  def("netCDisksSample", runAlgSample<anomaly::CenterNetDisks>);
  def("cDisksSample", runAlgSampleC);
  

}
