/* Copyright (C) University of Utah - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Michael Matheny <michaelmathen@gmail.com>, 2015-2016
 */
#ifndef ANOMALY_REGION_SCAN
#define ANOMALY_REGION_SCAN
#include <vector>
#include "Region.hpp"
#include "Point.hpp"
#include "Utility.hpp"

namespace anomaly {
  using namespace std;

  double norm(Point const& pt) {
    return pt.getX() * pt.getX() + pt.getY() * pt.getY(); 
  }

  struct AllDisks {
    template <typename RI>
    Disk operator()(RI nB, RI nE,
		    int netSize,
		    int sampleSize,
		    double rho) {
      vector<Point> netS;
      vector<Point> aS;
      vector<Point> bS;    
      dualSample(nB, nE, netS, aS, bS, netSize, sampleSize);
      return run(netS.begin(), netS.end(), aS.begin(), aS.end(), bS.begin(), bS.end(), rho);
    }
    
    template <typename RI>
    Disk run(RI nB, RI nE,
	     RI asBegin, RI asEnd,
	     RI bsBegin, RI bsEnd,
	     double rho) {

      //Remove any duplicates.
      sort(nB, nE, [](Point const& pt1, Point const& pt2) {
	  return norm(pt1) < norm(pt2);
	});
      auto nEnd = unique(nB, nE, [](Point const& pt1, Point const& pt2) {
	  return pt1.sameLoc(pt2);
	});
      auto nBegin = nB;
      
      vector<Point> netSampleSorted(nBegin, nEnd);

      Disk currMax;
      currMax.setNumPoints(0, 0);
      currMax.setNumAnomalies(0, 0);
      vector<int> aCountsA(nEnd - nBegin, 0);
      vector<int> bCountsA(nEnd - nBegin, 0);
      vector<int> aCountsR(nEnd - nBegin, 0);
      vector<int> bCountsR(nEnd - nBegin, 0);
      for (auto i = nBegin; i != nEnd - 1; i++) {
	for (auto j = i + 1; j != nEnd; j++) {
	  //Create a vector between the two points
	  double orthoX, orthoY;
	  findPerpVect(*i, *j, &orthoX, &orthoY);
	  double cX = (i->getX() + j->getX()) / 2.0;
	  double cY = (i->getY() + j->getY()) / 2.0;
	  
	  auto isNotCol = [&i, &j](Point const& pt) {
	    return !colinear(*i, *j, pt);
	  };
	  
	  // Partition these into a set of adding points and removing points
	  auto partitionF = [orthoX, orthoY, cX, cY](Point const& pt) {
	    return (pt.getX() - cX) * orthoX + (pt.getY() - cY) * orthoY <= 0;
	  };
	  
	  auto orderF = [orthoX, orthoY, &i, &j, cX, cY](Point const& pt) {
	    // If the point lines up with either of the reference
	    // point then we take this to be a disk defined by only
	    // the reference points.
	    if (colinear(*i, *j, pt)){
	      return numeric_limits<double>::infinity();
	    } else {
	      // We are projecting a vector created between
	      //the disk center and center point between the two points.
	      double a, b;
	      solveCircle3(*i, *j, pt, a, b);
	      return orthoX * (a - cX) + orthoY * (b - cY);
	    }
	    
	  };
	  
	  auto compF = [&orderF] (Point const& pt1, Point const& pt2) {
	    return orderF(pt1) < orderF(pt2);
	  };
	  auto onLine = [i, j](Point const& pt) {
	    return onLineSegment(*i, *j, pt);
	  };

	  auto asIterEnd = partition(asBegin, asEnd, isNotCol);
	  auto bsIterEnd = partition(bsBegin, bsEnd, isNotCol);
	  auto aCount = count_if(asIterEnd, asEnd, onLine);
	  auto bCount = count_if(bsIterEnd, bsEnd, onLine);
	  auto nIterEnd = partition(netSampleSorted.begin(), netSampleSorted.end(), isNotCol);
	  sort(netSampleSorted.begin(), nIterEnd, compF);
	  auto aHigherIt = partition(asBegin, asIterEnd, partitionF);
	  auto bHigherIt = partition(bsBegin, bsIterEnd, partitionF);

	  fill(aCountsR.begin(), aCountsR.end(), 0);
	  fill(bCountsR.begin(), bCountsR.end(), 0);
	  fill(aCountsA.begin(), aCountsA.end(), 0);
	  fill(bCountsA.begin(), bCountsA.end(), 0);
	  /*Probably most of the time is spent here*/
	  partial_counts(asBegin, aHigherIt,
			 netSampleSorted.begin(), nIterEnd,
			 aCountsR, compF);
	  aCount = aHigherIt - asBegin + aCount;
	  partial_counts(bsBegin, bHigherIt,
			 netSampleSorted.begin(), nIterEnd,
			 bCountsR, compF);
	  bCount = bHigherIt - bsBegin + bCount;
	  partial_counts(aHigherIt, asIterEnd,
			 netSampleSorted.begin(), nIterEnd,
			 aCountsA, compF);
	  partial_counts(bHigherIt, bsIterEnd,
			 netSampleSorted.begin(), nIterEnd,
			 bCountsA, compF);
	  /*----------------------------------------------*/
	  
	  auto size = nIterEnd - netSampleSorted.begin();
	  for (int k = 0; k < size; k++) {
	    aCount += aCountsA[k] - aCountsR[k];
	    bCount += bCountsA[k] - bCountsR[k];	    
	    double b_hat = bCount / static_cast<double>(bsEnd - bsBegin);
	    double a_hat = aCount / static_cast<double>(asEnd - asBegin);
	    
	    //This is not conservative at all
	    Disk currDisk(*i, *j, netSampleSorted[k]);
	    currDisk.setNumAnomalies(aCount, asEnd - asBegin);
	    currDisk.setNumPoints(bCount, bsEnd - bsBegin);
	    if (currMax.statistic() <= currDisk.statistic()) {
	      currMax = currDisk;
	    }
	  }
	}
      }
      return currMax;
    }
  };

  

  struct MoreDisks {
    template <typename RI>
    Disk operator()(RI nB, RI nE,
		    int netSize,
		    int sampleSize,
		    double rho) {
      vector<Point> netS;
      vector<Point> aS;
      vector<Point> bS;    
      dualSample(nB, nE, netS, aS, bS, netSize, sampleSize);
      return run(netS.begin(), netS.end(), aS.begin(), aS.end(), bS.begin(), bS.end(), rho);
    }

    template <typename T, typename F>
    void sort2(T& t, F& f) {
      if (f(**t[1]) < f(**t[0]))
	swap(t[0], t[1]);
    }
    
    template <typename T, typename F>
    void sort3(T& t, F& f) {
      sort2(t, f);
      if (f(**t[2]) < f(**t[0])) {
	swap(t[2], t[0]);
	swap(t[2], t[1]);
      } else if (f(**t[2]) < f(**t[1])) {
	swap(t[2], t[1]);
      }

    }
    
    template <typename T, typename F>
    void sort4(T& t, F& f) {
      sort2(t, f);
      auto l = t + 2;
      sort2(l, f);
      if (f(**t[2]) < f(**t[0])) {
	swap(t[2], t[0]);
	swap(t[1], t[2]);
	if (f(**t[3]) < f(**t[2])) {
	  swap(t[3], t[2]);
	}
      } else {
	if (f(**t[2]) < f(**t[1])) {
	  swap(t[2], t[1]);
	  if (f(**t[3]) < f(**t[2])) {
	    swap(t[3], t[2]);
	  }
	}
      }
    }

    
    template <typename RI>
    Disk run(RI nB, RI nE,
	     RI asBegin, RI asEnd,
	     RI bsBegin, RI bsEnd,
	     double rho) {
      
      //Remove any duplicates.
      sort(nB, nE, [](Point const& pt1, Point const& pt2) {
	  return norm(pt1) < norm(pt2);
	});
      auto nEnd = unique(nB, nE, [](Point const& pt1, Point const& pt2) {
	  return pt1.sameLoc(pt2);
	});
      auto nBegin = nB;
      
      Disk currMax;
      currMax.setNumPoints(0, 0);
      currMax.setNumAnomalies(0, 0);

      for (auto i = nBegin; i != nEnd - 1; i++) {
	for (auto j = i + 1; j != nEnd; j++) {
	  //Create a vector between the two points
	  double orthoX, orthoY;
	  findPerpVect(*i, *j, &orthoX, &orthoY);
	  double cX = (i->getX() + j->getX()) / 2.0;
	  double cY = (i->getY() + j->getY()) / 2.0;
	  
	  auto isNotCol = [&i, &j](Point const& pt) {
	    return !colinear(*i, *j, pt);
	  };
	  
	  // Partition these into a set of adding points and removing points
	  auto partitionF = [orthoX, orthoY, cX, cY](Point const& pt) {
	    return (pt.getX() - cX) * orthoX + (pt.getY() - cY) * orthoY <= 0;
	  };
	  
	  auto orderF = [orthoX, orthoY, &i, &j, cX, cY](Point const& pt) {
	    // If the point lines up with either of the reference
	    // point then we take this to be a disk defined by only
	    // the reference points.
	    if (colinear(*i, *j, pt)){
	      return numeric_limits<double>::infinity();
	    } else {
	      // We are projecting a vector created between
	      //the disk center and center point between the two points.
	      double a, b;
	      solveCircle3(*i, *j, pt, a, b);
	      return orthoX * (a - cX) + orthoY * (b - cY);
	    }
	  };
	  
	  auto compF = [&orderF] (Point const& pt1, Point const& pt2) {
	    return orderF(pt1) < orderF(pt2);
	  };
	  auto onLine = [i, j](Point const& pt) {
	    return onLineSegment(*i, *j, pt);
	  };

	  auto asIterEnd = partition(asBegin, asEnd, isNotCol);
	  auto bsIterEnd = partition(bsBegin, bsEnd, isNotCol);
	  auto aCount = count_if(asIterEnd, asEnd, onLine);
	  auto bCount = count_if(bsIterEnd, bsEnd, onLine);
	  
	  auto aHigherIt = partition(asBegin, asIterEnd, partitionF);
	  sort(asBegin, aHigherIt, compF);
	  sort(aHigherIt, asIterEnd, compF);
	  auto bHigherIt = partition(bsBegin, bsIterEnd, partitionF);
	  sort(bsBegin, bHigherIt, compF);
	  sort(bHigherIt, bsIterEnd, compF);
	  auto aRemov = asBegin, aAdd = aHigherIt;
	  auto bRemov = bsBegin, bAdd = bHigherIt;
	  while (true) {
	    
	    decltype(&aRemov) iters[4];
	    int eli = 0;
	    if (aRemov != aHigherIt) {
	      iters[eli] = &aRemov;
	      eli++;
	    }
	    if (bRemov != bHigherIt) {
	      iters[eli] = &bRemov;
	      eli++;
	    }
	    if (aAdd != asEnd) {
	      iters[eli] = &aAdd;
	      eli++;
	    }
	    if (bAdd != bsEnd) {
	      iters[eli] = &bAdd;
	      eli++;
	    }

	    decltype(&aRemov) k;
	    if (eli == 0) {
	      break;
	    } else if (eli == 1) {
	      k = iters[0];
	    } else {
	      if (eli == 2) {
		sort2(iters, orderF);
	      } else if (eli == 3) {
		sort3(iters, orderF);
	      } else 
		sort4(iters, orderF);
	      k = iters[0];
	      if (orderF(**k) == orderF(**iters[1])) {
		(*k)++;
		continue;
	      }
	    }
	    double b_hat = static_cast<double>(bAdd - bRemov + bCount) / (bsEnd - bsBegin);
	    double a_hat = static_cast<double>(aAdd - aRemov + aCount) / (asEnd - asBegin);
	    double eps = 2 / static_cast<double>(nEnd - nBegin);
	    
	    if (rhoInRange(a_hat, b_hat, rho, eps)) {
	      Disk currDisk(*i, *j, **k);
	      currDisk.setNumAnomalies(aAdd - aRemov + aCount, asEnd - asBegin);
	      currDisk.setNumPoints(bAdd - bRemov + bCount, bsEnd - bsBegin);
	      if (currMax.statistic() <= currDisk.statistic()) {
		currMax = currDisk;
	      }
	    }
	    ++(*k);
	  }
	}
      }
      return currMax;
    }
  };


  /////////////////////////
  //End of All Disks//
  /////////////////////////
  

  struct CenteredDisks {

    template <typename RI>
    Disk operator()(RI nB, RI nE,
		    int sampleSize,
		    double rho) {
      vector<Point> aS;
      vector<Point> bS;    
      singleSample(nB, nE, aS, bS, sampleSize);
      return run(aS.begin(), aS.end(), bS.begin(), bS.end(), rho);
    }
    
    template <typename RI>
    Disk run(RI aSB, RI aSE,
	     RI bSB, RI bSE,
	     double rho) {

      vector<Point> centers;

      centers.insert(centers.end(), aSB, aSE);
      centers.insert(centers.end(), bSB, bSE);

      Disk currMax;
      currMax.setNumPoints(0, centers.size());
      currMax.setNumAnomalies(0, centers.size());

      for (auto const& pt : centers) {
	auto compF = [&pt](const Point& pt1, const Point& pt2) {
	       return pt.distanceSquared(pt1) < pt.distanceSquared(pt2);
	};
	sort(aSB, aSE, compF);
	sort(bSB, bSE, compF);
	auto ia = aSB;
	auto ib = bSB;
	bool flag = true;
	do {
	  Disk currDisk;
	  if (ia == aSE) {
	    ib = bSE;
	    currDisk = Disk(pt, *(bSE - 1));
	    flag = false;
	  } else if (ib == bSE) {
	    ia = aSE;
	    currDisk = Disk(pt, *(aSE - 1));
	    flag = false;
	  } else if (pt.distanceSquared(*ia) < pt.distanceSquared(*ib)) {
	    currDisk = Disk(pt, *ia);
	    ia++;
	  } else {
	    currDisk = Disk(pt, *ib);
	    ib++;
	  }
	  double b_hat = (ib - bSB) / static_cast<double>(bSE - bSB);
	  double a_hat = (ia - aSB) / static_cast<double>(aSE - aSB);
	  //This is not conservative at all
	  double eps = 2 / static_cast<double>(bSE - bSB);
	  if (rhoInRange(a_hat, b_hat, rho, eps)) {
	    currDisk.setNumAnomalies(ia - aSB, aSE - aSB);
	    currDisk.setNumPoints(ib - bSB, bSE - bSB);
	    if (currMax.statistic() <= currDisk.statistic()) {
	      currMax = currDisk;
	    }
	  }
	} while (flag);
      }
      return currMax;
    }
  };

  struct CenterNetDisks {
    template<typename RI>
    Disk operator()(RI nB, RI nE,
		    int netSize,
		    int sampleSize,
		    double rho) {
      vector<Point> netS;
      vector<Point> aS;
      vector<Point> bS;
      dualSample(nB, nE, netS, aS, bS, netSize, sampleSize);
      return run(netS.begin(), netS.end(), aS.begin(), aS.end(), bS.begin(), bS.end(), rho);
    }
    
    template <typename RI>
    Disk run(RI nB, RI nE,
	     RI aSB, RI aSE,
	     RI bSB, RI bSE,
	     double rho) {
      
      Disk currMax;
      currMax.setNumPoints(0, bSE - bSB);
      currMax.setNumAnomalies(0, aSE - aSB);

      for (auto cp = nB; cp != nE; ++cp) {
	auto const& pt = *cp;
	auto compF = [&pt](const Point& pt1, const Point& pt2) {
	  return pt.distanceSquared(pt1) < pt.distanceSquared(pt2);
	};
	sort(aSB, aSE, compF);
	sort(bSB, bSE, compF);
	auto ia = aSB;
	auto ib = bSB;
	bool flag = true;
	do {
	  Disk currDisk;
	  if (ia == aSE) {
	    ib = bSE;
	    currDisk = Disk(pt, *(bSE - 1));
	    flag = false;
	  } else if (ib == bSE) {
	    ia = aSE;
	    currDisk = Disk(pt, *(aSE - 1));
	    flag = false;
	  } else if (pt.distanceSquared(*ia) < pt.distanceSquared(*ib)) {
	    currDisk = Disk(pt, *ia);
	    ia++;
	  } else {
	    currDisk = Disk(pt, *ib);
	    ib++;
	  }
	  double b_hat = (ib - bSB) / static_cast<double>(bSE - bSB);
	  double a_hat = (ia - aSB) / static_cast<double>(aSE - aSB);
	  //This is not conservative at all
	  double eps = 2 / static_cast<double>(nE - nB);
	  if (rhoInRange(a_hat, b_hat, rho, eps)) {
	    currDisk.setNumAnomalies(ia - aSB, aSE - aSB);
	    currDisk.setNumPoints(ib - bSB, bSE - bSB);
	    if (currMax.statistic() <= currDisk.statistic()) {
	      currMax = currDisk;
	    }
	  }
	} while (flag);
      }
      return currMax;
    }
  };

  struct SatScan {
    template <typename RI>
    Disk operator()(RI nB, RI nE) {
      return run(nB, nE);
    }

    template <typename RI>
    Disk run(RI nB, RI nE) {
      Disk currMax;
      currMax.setNumPoints(0, 0);
      currMax.setNumAnomalies(0, 0);
      vector<Point> sortedPoints(nB, nE);
      auto totalMeasuredCount = count_if(nB, nE, [](Point const& pt) {
	  return pt.getLabel();
	});
      for (auto i = nB; i != nE - 1; i++) {
	for (auto j = i + 1; j != nE; j++) {
	  //Create a vector between the two points
	  double orthoX, orthoY;
	  findPerpVect(*i, *j, &orthoX, &orthoY);
	  double cX = (i->getX() + j->getX()) / 2.0;
	  double cY = (i->getY() + j->getY()) / 2.0;
	  
	  auto isNotCol = [&i, &j](Point const& pt) {
	    return !colinear(*i, *j, pt);
	  };
	  
	  // Partition these into a set of adding points and removing points
	  auto partitionF = [orthoX, orthoY, cX, cY](Point const& pt) {
	    return (pt.getX() - cX) * orthoX + (pt.getY() - cY) * orthoY <= 0;
	  };
	  
	  auto orderF = [orthoX, orthoY, &i, &j, cX, cY](Point const& pt) {
	    // If the point lines up with either of the reference
	    // point then we take this to be a disk defined by only
	    // the reference points.
	    if (colinear(*i, *j, pt)){
	      return numeric_limits<double>::infinity();
	    } else {
	      // We are projecting a vector created between
	      //the disk center and center point between the two points.
	      double a, b;
	      solveCircle3(*i, *j, pt, a, b);
	      return orthoX * (a - cX) + orthoY * (b - cY);
	    }
	  };
	  
	  auto compF = [&orderF] (Point const& pt1, Point const& pt2) {
	    return orderF(pt1) < orderF(pt2);
	  };
	  auto onLine = [i, j](Point const& pt) {
	    return onLineSegment(*i, *j, pt);
	  };
	  
	  auto iterEnd = partition(sortedPoints.begin(), sortedPoints.end(), isNotCol);

	  auto measuredOnLine = count_if(iterEnd, sortedPoints.end(), [&onLine](Point const& pt) {
	      return pt.getLabel() && onLine(pt);
	    });
	  auto baseOnLine = count_if(iterEnd, sortedPoints.end(), [&onLine](Point const& pt) {
	      return onLine(pt);
	    });

	  auto higherIt = partition(sortedPoints.begin(), iterEnd, partitionF);
	  sort(sortedPoints.begin(), higherIt, compF);
	  sort(higherIt, iterEnd, compF);
	  
	  auto aCount = count_if(sortedPoints.begin(), higherIt, [](Point const& pt) {
	      return pt.getLabel();
	    });
	  auto bCount = higherIt - sortedPoints.begin();
	  auto remov = sortedPoints.begin(), add = higherIt;
	  
	  while (true) {
	    decltype(i) k;
	    if (add == iterEnd && remov == higherIt) {
	      break;
	    } else if (add == iterEnd) {
	      k = remov;
	    } else if (remov == higherIt) {
	      k = add;
	    } else if (orderF(*add) < orderF(*remov)) {
	      k = add;
	    } else {
	      k = remov;
	    }
	    Disk currDisk(*i, *j, *k);
	    currDisk.setNumAnomalies(aCount + measuredOnLine, totalMeasuredCount);
	    currDisk.setNumPoints(bCount + baseOnLine, sortedPoints.size());
	    if (currMax.statistic() <= currDisk.statistic()) {
	      currMax = currDisk;
	    }
	    
	    if (add == iterEnd) {
	      aCount -= add->getLabel();
	      bCount--;
	      remov++;
	    } else if (remov == higherIt) {
	      aCount += add->getLabel();
	      bCount++;
	      add++;
	    } else if (orderF(*add) < orderF(*remov)) {
	      aCount += add->getLabel();
	      bCount++;
	      add++;
	    } else {
	      aCount -= add->getLabel();
	      bCount--;
	      remov++;
	    }
	  }
	}
      }
      return currMax;
    }
  };

}

#endif
