#ifndef RECT_SCAN
#define RECT_SCAN

#include <map>
#include "Region.hpp"

#include "Utility.hpp"
#include "Point.hpp"
//#include "OrthoCount.hpp"

namespace anomaly {
  struct RectangleScan {

    template <typename T, typename Struct>
    auto less_or_equal(T key, Struct& el) -> decltype(el.begin()) {
      auto lb = el.lower_bound(key);
      return lb->first == key ? lb : --lb;
    }

    /*
      We assume that all containers are sorted by y axis.
     */
    template <typename RI>
    vector<vector<double>> generateStrips(RI begin, RI end,
					  RI sBegin, RI sEnd) {
      vector<vector<double>> strips(end - begin - 1, vector<double>());
      auto curr_strip = strips.begin();
      auto sb = sBegin;
      //Find the begining of the first point
      for (;sb != sEnd && sb->getY() < begin->getY(); sb++);
      //Do all the elements between.
      for (auto i = begin + 1; i != end; ++i) {
	for (; sb != sEnd && sb->getY() < i->getY(); ++sb) {
	  curr_strip->push_back(sb->getX());
	}
	sort(curr_strip->begin(), curr_strip->end());
	++curr_strip;
      }
      return strips;
    }

    void insertNode(Point const& pt, map<double, vector<double>>& nodes) {
      //Splits the preceading node in the heap from the current node we are inserting.
      auto lb = nodes.lower_bound(pt.getX());
      if (lb == nodes.end() || lb->first != pt.getX()) {
	lb--;
	auto unaryF = [&pt](double el) {
	  return el < pt.getX();
	};
	auto elsAfterPart = partition(lb->second.begin(), lb->second.end(), unaryF);
	if (lb->second.end() - elsAfterPart > 0) {
	  nodes.insert(make_pair(pt.getX(), vector<double>(elsAfterPart, lb->second.end())));
	  lb->second.erase(elsAfterPart);
       } else {
	  nodes.insert(make_pair(pt.getX(), vector<double>()));
       }
      }
    }

    template <typename T>
    void updateNodes(T b, T e, map<double, vector<double>>& nodes) {
      // We run through all the nodes before the last one.
      auto nI = ++nodes.begin();
      auto nIB = nodes.begin();
      for (; nI != nodes.end(); ++nI) {
	for (; b != e && *b < nI->first; ++b) {
	  nIB->second.push_back(*b);
	}
	nIB++;
      }
      // Now we finish with the last node
      for (; b != e; ++b) {
	nIB->second.push_back(*b);
      }
    }

    template<typename T>
    int sumBox(T l_it, T r_it) {
      int count = 0;
      for (;l_it != r_it; ++l_it) {
	count += l_it->second.size();
      }
      return count;
    }
    
    template <typename RI>
    Rectangle operator()(RI nB, RI nE,
			 int netSize,
			 int sampleSize,
			 double rho) {
      vector<Point> netS;
      vector<Point> aS;
      vector<Point> bS;    
      dualSample(nB, nE, netS, aS, bS, netSize, sampleSize);
      return run(netS.begin(), netS.end(), aS.begin(), aS.end(), bS.begin(), bS.end(), rho);
    }


    /*
      This should run in roughly O(slogs + ns + n^4) time.
     */
    using node_t = tuple<vector<Point*>, vector<Point*>>;
    template <typename RI>
    Rectangle run(RI begin, RI end,
		  RI asBegin, RI asEnd,
		  RI bsBegin, RI bsEnd,
		  double rho) {
      auto cmpX = [](Point const& pt1, Point const& pt2) {
	return pt1.getX() < pt2.getX();
      };
      
      auto cmpY = [](Point const& pt1, Point const& pt2) {
	return pt1.getY() < pt2.getY();
      };
      sort(begin, end, cmpY);
      sort(asBegin, asEnd, cmpY);
      sort(bsBegin, bsEnd, cmpY);

      //At the end of this part we have each 2n strips where each strip is sorted by x axis
      auto aStrips = generateStrips(begin, end, asBegin, asEnd);
      auto bStrips = generateStrips(begin, end, bsBegin, bsEnd);

      Rectangle maxRect;
      maxRect.setNumAnomalies(0, 0);
      maxRect.setNumPoints(0, 0);
      double maxScan = maxRect.statistic();
      auto ai = aStrips.begin();
      auto bi = bStrips.begin();
      for (auto i = begin;i != end - 1; ai++, bi++, i++) {
	map<double, vector<double>> aNodes;
	map<double, vector<double>> bNodes;
	aNodes.insert(make_pair(-numeric_limits<double>::infinity(), vector<double>()));
	bNodes.insert(make_pair(-numeric_limits<double>::infinity(), vector<double>()));
	insertNode(*i, aNodes);
	insertNode(*i, bNodes);
	auto aj = ai;
	auto bj = bi;
	for (auto j = i + 1; j != end; aj++, bj++, j++) {
	  insertNode(*j, aNodes);
	  insertNode(*j, bNodes);
	  updateNodes(aj->begin(), aj->end(), aNodes);
	  updateNodes(bj->begin(), bj->end(), bNodes);

	  double left_pt = i->getX() < j->getX() ? i->getX() : j->getX();
	  double right_pt = i->getX() >= j->getX() ? i->getX() : j->getX();

	  auto l_itA = aNodes.lower_bound(left_pt);
	  auto r_itA = aNodes.lower_bound(right_pt);
	  auto l_itB = bNodes.lower_bound(left_pt);
	  auto r_itB = bNodes.lower_bound(right_pt);
	  int aCount = sumBox(l_itA, r_itA);
	  int bCount = sumBox(l_itB, r_itB);
	  #ifdef DEBUG
	  cout << "bottom left = " << left_pt << " " << i->getY() << endl;
	  cout << "top right side = " << right_pt << " " << j->getY() << endl;
	  cout << "aCount = " << aCount << endl;
	  cout << "bCount = " << bCount << endl;
	  #endif

	  // Now we can iterate over every rectangle that has *i and *j as its upper and lower
	  // elements respectively
	  // We could do every pair, but that would overcount the rectangles so we actually want
	  // to do every pair that lies outside of i->getX() and j->getX().
	  // ---------------*i--------------------------
	  //                      3         
	  // 1        2
	  //                                     4
	  //-----------------------------*j-------------
	  //So above we do rectangles (1, 4, *i, *j) and (2, 4, *i, *j), but not (2, 3, *i, *j)
	  //Now we compute all rectangles defined by 3 and 4 points.

	  do {
	    int aCount_tmp = aCount;
	    int bCount_tmp = bCount;
	    auto ra = r_itA;
	    auto rb = r_itB;
	    do {
	      double b_hat = bCount_tmp / static_cast<double>(bsEnd - bsBegin);
	      double a_hat = aCount_tmp / static_cast<double>(asEnd - asBegin);
	      //This is not conservative at all
	      double eps = 2 / static_cast<double>(end - begin);
	      if (rhoInRange(a_hat, b_hat, rho, eps)) {
		Rectangle rect(l_itA->first, ra->first, i->getY(), j->getY());
		rect.setNumAnomalies(aCount_tmp, asEnd - asBegin);
		rect.setNumPoints(bCount_tmp, bsEnd - bsBegin);
		double newSt = rect.statistic();
		if (newSt >= maxScan) {
		  maxRect = rect;
		  maxScan = newSt;
		}
	      }
	      aCount_tmp += ra->second.size();
	      bCount_tmp += rb->second.size();
	      ra++;
	      rb++;
	    } while (ra != aNodes.end());
	    l_itA--;
	    l_itB--;
	    aCount += l_itA->second.size();
	    bCount += l_itB->second.size();
	  } while (l_itA != aNodes.begin());
	} 
      }
      return maxRect;
    }
  };
}

#endif
