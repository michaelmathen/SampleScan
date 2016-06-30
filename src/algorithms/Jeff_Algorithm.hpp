/* Copyright (C) University of Utah - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Michael Matheny <michaelmathen@gmail.com>, 2015-2016
 */
#ifndef JEFF_ALGORITHM_MAIN
#define JEFF_ALGORITHM_MAIN

#include <vector>
#include <time.h>
#include <memory>

#include "Point.hpp"

extern "C" {
   #include "algorithms/old/dataset.h"
   #include "algorithms/old/function.h"
}

namespace anomaly{

  using namespace std;
  /* The function to imitate Jeff's code
   * Para:
   * 	points, the points array
   * 	output, the file name of output file
   * 	eps, the epsilon value
   */


  struct JeffAlgo {
    template<typename RI>
    Rectangle operator() (RI begin, RI end, double eps) {
      double x_min, x_max, y_min, y_max, maxnr, maxnb;
      ad_fxn A;
      data_set ds;
      shared_ptr<data_point> dsData;
      d_fxn D;

      ds.npts = end - begin;
      dsData = shared_ptr<data_point>(static_cast<data_point*>(malloc(ds.npts * sizeof(data_point))), free);
      ds.data = dsData.get();
      ds.bpts = ds.rpts = 0;
      for (int i = 0; i < ds.npts; ++i){
      	ds.data[i].x = begin[i].getX();
      	ds.data[i].y = begin[i].getY();
      	ds.data[i].blue = 1.0;
      	ds.bpts += ds.data[i].blue;
      	ds.data[i].red = begin[i].getLabel() ? 1.0 : 0.0;
      	ds.rpts += ds.data[i].red;
      }
      // find approximate function
      set_Poisson_fxn(&D);
      smart_gridding_fxn(&D, &A, eps, ds.npts);
      create_heap_ds(&ds);

      Rectangle currMax;
      for (int i = 0; i < A.n; i++) {
	set_red_blue_weight_ds(&ds, A.x[i]/ds.rpts, A.y[i]/ds.bpts);
	double m = max_discrepancy_ds(&ds, &x_min, &x_max, &y_min, &y_max, &maxnr, &maxnb);
	Rectangle currRegion(x_min, x_max, y_min, y_max);
	currRegion.useKStat(m);
	if (currRegion.statistic() >= currMax.statistic()) {
	  currMax = currRegion;
	}
      }
      return currMax;
    }
  };
}

#endif
