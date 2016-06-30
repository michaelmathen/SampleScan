/* Copyright (C) University of Utah - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Michael Matheny <michaelmathen@gmail.com>, 2015-2016
 */
#ifndef ANOMALY_UTILITIES
#define ANOMALY_UTILITIES
#include <cmath>
#include <cstddef>
#include <vector>
#include <random>
#include <fstream>
#include "Point.hpp"

#include <boost/algorithm/string.hpp>
namespace anomaly {


  template <typename T, typename C>
  inline void partial_counts(T begin, T end,
			     T break_begin, T break_end, // Assumed to be sorted
			     vector<int>& counts,
			     C compF) {
    //Partitions based on the break points.
    for (; begin != end; begin++) {
      auto lb = lower_bound(break_begin, break_end, *begin, compF);
      counts[lb - break_begin] += 1;
    }
  }
  
  inline bool rhoInRange(double m, double b, double rho, double eps) {
    double alpha = exp(- 1 / rho);
    return (m >= alpha + eps) && (m <= 1 - alpha - eps) &&
      (b >= rho + eps) && (b <= 1 - rho - eps);
  }

  // struct CompF {
    
  //   double a0, b0, a, b;
  //   CompF(Point const& pt1, Point const& pt2) {
  //     double x1, x2, y1, y2;
  //     pt1.getLoc(x1, y1);
  //     pt2.getLoc(x2, y2);
  //     a = y2 - y1;
  //     b = x1 - x2;
  //     b0 = (y1 + y2) / 2;
  //     a0 = (x1 + x2) / 2;
  //   }
  //   inline bool operator()(Point const& pt1, Point const& pt2) const {
  //     double x1, x2, y1, y2;
  //     pt1.getLoc(x1, y1);
  //     pt2.getLoc(x2, y2);
  //     a = y2 - y1;
  //     b = x1 - x2;
  //     b0 = (y1 + y2) / 2;
  //     a0 = (x1 + x2) / 2;
      
  //   }
  // };

  /*
    Takes three points that are assumed to 
    lie on the boundary of a disk.
    INPUTS
    pt1, pt2, pt3 -- Three points on the boundary of the disk.
    (a, b)        -- Center of the defined disk.
    r             -- Radius of the defined disk.
   */
  void solveCircle3(Point const& pt1,
		    Point const& pt2,
		    Point const& pt3,
		    double &a,
		    double &b,
		    double &r);

  bool colinear(Point const& pt1,
		Point const& pt2,
		Point const& pt3); 
    /*
    Takes three points that are assumed to 
    lie on the boundary of a disk.
    INPUTS
    pt1, pt2, pt3 -- Three points on the boundary of the disk.
    (a, b)        -- Center of the defined disk, but doesn't scale by 
                     the determinant

   */
  void solveCircle3(Point const& pt1,
			   Point const& pt2,
			   Point const& pt3,
			   double &a,
			   double &b);

  /*
    Checks to see if pt3 is on the line segment between pt1 and pt2
  */
  bool onLineSegment(Point const& pt1,
		     Point const& pt2,
		     Point const& pt3);
  /*
    Finds the center and radius of a circle defined by two oposing 
    points.
    INPUTS
    pt1, pt2 -- Opposiing points
    RETURNS
    (a, b) -- Center of the disk
    r      -- The radius of the disk
  */
  void solveCircle2(Point const& pt1,
		    Point const& pt2,
		    double &a,
		    double &b,
		    double &r);

  bool inDisk(double a, double b, double r, Point const& pt);

  /*
    Taken from google code
   */
  /*=====================================================================*
   *                   Copyright (C) 2011 Paul Mineiro                   *
   * All rights reserved.                                                *
   *                                                                     *
   * Redistribution and use in source and binary forms, with             *
   * or without modification, are permitted provided that the            *
   * following conditions are met:                                       *
   *                                                                     *
   *     * Redistributions of source code must retain the                *
   *     above copyright notice, this list of conditions and             *
   *     the following disclaimer.                                       *
   *                                                                     *
   *     * Redistributions in binary form must reproduce the             *
   *     above copyright notice, this list of conditions and             *
   *     the following disclaimer in the documentation and/or            *
   *     other materials provided with the distribution.                 *
   *                                                                     *
   *     * Neither the name of Paul Mineiro nor the names                *
   *     of other contributors may be used to endorse or promote         *
   *     products derived from this software without specific            *
   *     prior written permission.                                       *
   *                                                                     *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND              *
   * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,         *
   * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES               *
   * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE             *
   * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER               *
   * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,                 *
   * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES            *
   * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE           *
   * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR                *
   * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF          *
   * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT           *
   * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY              *
   * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE             *
   * POSSIBILITY OF SUCH DAMAGE.                                         *
   *                                                                     *
   * Contact: Paul Mineiro <paul@mineiro.com>                            *
   *=====================================================================*/
  inline float fastlog2 (float x) {
    union { float f; uint32_t i; } vx = { x };
    union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
    float y = vx.i;
    y *= 1.1920928955078125e-7f;

    return y - 124.22551499f
      - 1.498030302f * mx.f 
      - 1.72587999f / (0.3520887068f + mx.f);
  }

  inline float fast_log (float x) {
    return 0.69314718f * fastlog2 (x);
  }
  
  inline float fkulldorff(int cIn, int cTotal, int aIn, int aTotal) {
    /*
      Calculate the kulldorff scan statistic. 
    */
    if (cIn == 0 || cTotal == 0 ||
	cIn == cTotal || aTotal == 0) {
      return 0.0;
    }
    float br = (float)cIn / cTotal;
    if (aIn == 0) {
      return log(1 / (1 - br));
    }
    if (aIn == aTotal) {
      return log(1 / br);
    }
  
    float mr = (float)aIn / aTotal;
    return mr * log(mr / br) + (1 - mr) * log((1 - mr) / (1 - br));
  }

  
  inline double kulldorff(int cIn, int cTotal, int aIn, int aTotal) {
    /*
      Calculate the kulldorff scan statistic. 
    */
    if (cIn == 0 || cTotal == 0 ||
	cIn == cTotal || aTotal == 0) {
      return 0.0;
    }
    double br = (double)cIn / cTotal;
    if (aIn == 0) {
      return log(1 / (1 - br));
    }
    if (aIn == aTotal) {
      return log(1 / br);
    }
  
    double mr = (double)aIn / aTotal;
    return mr * log(mr / br) + (1 - mr) * log((1 - mr) / (1 - br));
  }

  void findPerpVect(Point const& p1, Point const& p2, double* u, double* v);

  void findOrthoVector(double at1, double at2,
		       double &x1, double &y1);
  
  /*
    Uses reservior sampling to generate samples based on some function f
   */
  template <typename Value, typename Size, typename F>
  class StreamSample {
    vector<Value>& sample;
    Size sampleSize;
    Size __i;
    mt19937 generator;

    F func;
  public:
    StreamSample(Size sampleSize,
		 vector<Value>& vect,
		 F func) :
      sampleSize(sampleSize),
      __i(sampleSize),
      sample(vect),
      func(func)
    {
      sample.reserve(sampleSize);
      random_device rd;
      generator = mt19937(rd());
    }

    void update(Value const& currEl) {
      if (func(currEl)) {
	if (sample.size() < sampleSize) {
	  sample.push_back(currEl);
	}
	else {
	  uniform_int_distribution<ptrdiff_t> distribution(0, __i - 1);
	  auto ix = distribution(generator);
	  if (ix < sampleSize) {
	    sample[ix] = currEl;
	  }
	  __i++;
	}
      }
    }
  };

  struct {
    bool operator()(Point const& p) { return true; }
  } fBaseline;
  struct {
    bool operator() (Point const& p) { return p.getLabel(); }
  } fAnomaly;

  template<typename InputIt>
  void dualSample(InputIt begin, InputIt end,
		  vector<Point>& netSample, 
		  vector<Point>& aSample, 
		  vector<Point>& bSample,
		  long netSize,
		  long sampleSize) {
    
    vector<Point> aNet;
    vector<Point> bNet;
    StreamSample<Point, long, decltype(fBaseline)> bN(netSize, bNet, fBaseline);
    StreamSample<Point, long, decltype(fAnomaly)> aN(netSize, aNet, fAnomaly);
    StreamSample<Point, long, decltype(fBaseline)> bS(sampleSize, aSample, fBaseline);
    StreamSample<Point, long, decltype(fAnomaly)> aS(sampleSize, bSample, fAnomaly);
    for (; begin != end; begin++) {
      auto el = *begin;
      aN.update(el);
      bN.update(el);
      aS.update(el);
      bS.update(el);
    }
    netSample.insert(netSample.end(), aNet.begin(), aNet.end());
    netSample.insert(netSample.end(), bNet.begin(), bNet.end());
  }

  template<typename InputIt>
  void singleSample(InputIt begin, InputIt end,
		    vector<Point>& aSample, 
		    vector<Point>& bSample, 
		    long sampleSize) {
    StreamSample<Point, long, decltype(fBaseline)> bS(sampleSize, aSample, fBaseline);
    StreamSample<Point, long, decltype(fAnomaly)> aS(sampleSize, bSample, fAnomaly);
    for (; begin != end; begin++) {
      auto el = *begin;
      aS.update(el);
      bS.update(el);
    }
  }


  void parsePtFile(string const& fname, vector<Point>& pts);

  template<typename T>
  ostream& operator<<(ostream& os, vector<T>& v) {
    os << "[" ;
    for (auto val = v.begin(); val != v.end(); val++) {
      os << *val << ", ";
    }
    os << "]";
    return os;
  }


}
#endif
