#ifndef RECT_SCAN
#define RECT_SCAN
#include <omp.h>
#include <map>
#include "Region.hpp"

#include "Utility.hpp"
#include "Point.hpp"
//#include "OrthoCount.hpp"

namespace anomaly {

  template<typename T>
  struct Grid {
    using size_type = typename vector<T>::size_type;
    vector<T> grid;
    size_type n;
    Grid(size_type n) : grid(n * n, 0), n(n) {
      grid.resize(n * n);
    }
    T& operator()(size_type i, size_type j) {
      return grid[i * n + j];
    }
    T const& operator()(size_type i, size_type j) const {
      return grid[i * n + j];
    }
    size_type size() const {
      return n;
    }
  };

  template<typename P, typename RI>
  inline void fillGrid(P bx, P ex, P by, P ey, Grid<int>& counts, RI begin, RI end) {
    auto cmpX = [](Point const& pt1, Point const& pt2) {
      return pt1.getX() < pt2.getX();
    };
    auto cmpY = [](Point const& pt1, Point const& pt2) {
      return pt1.getY() < pt2.getY();
    };
    
    sort(begin, end, cmpX);
    // Divide the points in the x direction into strips
    // then sort each strip in the y direction.
    // Could be slightly sped up by splitting instead of sorting as this
    // will change from s * log(s) + s/n log(s/n) to s * log(n) + s/n * log(n)
    for (int i = 0; bx != ex; ++bx, ++i) { // for each strip
      auto e = find_if(begin, end, [&bx](Point const& el) {
	  return *bx < el.getX(); // Get the first element larger than this.
	});
      sort(begin, e, cmpY);
      int j = 0;
      for (auto tmpy = by; tmpy != ey; ++tmpy, ++j) { // for each grid cell
	while (begin->getY() <= *tmpy && begin != e) { // For each point.
	  counts(i, j) += 1;
	  begin++;
	}
      }
    }
    // Compute cummalative sums in the x direction
    for (int i = 1; i < counts.size(); i++) {
      for (int j = 0; j < counts.size(); j++) {
	counts(i, j) += counts(i - 1, j);	
      }
    }
    // Compute cummalative sums in the y direction
    for (int i = 0; i < counts.size(); i++) {    
      for (int j = 1; j < counts.size(); j++) {    
	counts(i, j) += counts(i, j - 1);
      }
    }

  }

  /*
    Calculates the grid cell by subtracting the appropriate cummulative sums
  */
  template<typename T>
  T computeCell(Grid<T> const& g, int ux, int uy, int lx, int ly) {
    return g(ux, uy) + g(lx, ly) - g(lx, uy) - g(ux, ly);
  }

  struct RectangleScanG {
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

    bool smaller(Rectangle r1, Rectangle r2) {
      return (r1.getUX() - r1.getLX()) * (r1.getUY() - r1.getLY()) > (r2.getUX() - r2.getLX()) * (r2.getUY() - r2.getLY());
    }
    /*
      This should run in roughly O(s logs + n^4) time.
     */
    template <typename RI>
    Rectangle run(RI begin, RI end,
		  RI asBegin, RI asEnd,
		  RI bsBegin, RI bsEnd,
		  double rho) {

      auto g = end - begin;
      vector<double> x_pos(g);
      vector<double> y_pos(g);      
      Grid<int> anomalies(g);
      Grid<int> baseline(g);
      int it = 0;
      for (auto b = begin; b != end; b++) {
	x_pos[it] = b->getX();
	y_pos[it] = b->getY();
	it++;
      }
      sort(x_pos.begin(), x_pos.end());
      sort(y_pos.begin(), y_pos.end());

      fillGrid(x_pos.begin(), x_pos.end(), y_pos.begin(), y_pos.end(),
	       anomalies, asBegin, asEnd);
      fillGrid(x_pos.begin(), x_pos.end(), y_pos.begin(), y_pos.end(),
	       baseline, bsBegin, bsEnd);

      const int aTotal = asEnd - asBegin;
      const int bTotal = bsEnd - bsBegin;
      const double eps = 2.0 / static_cast<double>(end - begin);
      const int num_threads = omp_get_max_threads();
      Rectangle rects[num_threads];
      #pragma omp parallel
      {
	int id = omp_get_thread_num();
	Rectangle maxRect;
	maxRect.setNumAnomalies(0, 0);
	maxRect.setNumPoints(0, 0);
	double maxScan = 0;
        #pragma omp for 
	for (int i = 0; i < g - 1; i++) { // left
	  for (int j = 0; j < g - 1; j++) { // top
	    for (int k = i + 1; k < g; k++) { // right
	      for (int m = j + 1; m < g; m++) { // bottom
		int aCount = computeCell(anomalies, i, j, k, m);
		int bCount = computeCell(baseline, i, j, k, m);
		double mr = (double) aCount / aTotal;
		double br = (double) bCount / bTotal;
		if (rhoInRange(mr, br, rho, eps)) {
		  Rectangle rect(x_pos[i], x_pos[k], y_pos[j], y_pos[m]);
		  rect.setNumAnomalies(aCount, aTotal);
		  rect.setNumPoints(bCount, bTotal);
		  double newSt = rect.statistic();
		  if (newSt > maxScan) {
		    maxRect = rect;
		    maxScan = newSt;
		  } else if (newSt == maxScan && smaller(maxRect, rect)) {
		    maxRect = rect;
		    maxScan = newSt;
		  }
		}
	      }
	    }
	  }
	}
	//return maxRect;
	rects[id] = maxRect;
      }
      Rectangle rect = rects[0];
      for (int i = 1; i < num_threads; i++) {
       	if (rects[i].statistic() > rect.statistic()) {
       	  rect = rects[i];
       	}
      }
      return rect;
    }
  };

  /*
    Try to figure out the maximum possible value that could occur if we continue.
   */
  inline double maxSub(double mr, double br, double eps, double rho) {
    rho = rho + eps;
    double alpha = exp(- 1 / rho) + eps;
    mr = min(1 - alpha, mr);
    br = min(1 - rho, br);
    mr = max(mr, alpha);
    br = max(br, rho);

    //Largest point occurs in one of the corners of (alpha, br), (mr, br),
    //(alpha, rho), (mr, rho)
    return max({alpha * log(alpha / br) + (1 - alpha) * log((1 - alpha) / (1 - br)),
	  alpha * log(alpha / rho) + (1 - alpha) * log((1 - alpha) / (1 - rho)),
	  mr * log(mr / rho) + (1 - rho) * log((1 - mr) / (1 - rho)),
	  mr * log(mr / br) + (1 - mr) * log((1 - mr) / (1 - br))});
  }
  
  struct RectangleScanGH {
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
    template <typename RI>
    Rectangle run(RI begin, RI end,
		  RI asBegin, RI asEnd,
		  RI bsBegin, RI bsEnd,
		  double rho) {

      auto g = end - begin;
      vector<double> x_pos(g);
      vector<double> y_pos(g);      
      Grid<int> anomalies(g);
      Grid<int> baseline(g);
      int it = 0;
      for (auto b = begin; b != end; b++) {
	x_pos[it] = b->getX();
	y_pos[it] = b->getY();
	it++;
      }
      sort(x_pos.begin(), x_pos.end());
      sort(y_pos.begin(), y_pos.end());

      fillGrid(x_pos.begin(), x_pos.end(), y_pos.begin(), y_pos.end(),
	       anomalies, asBegin, asEnd);
      fillGrid(x_pos.begin(), x_pos.end(), y_pos.begin(), y_pos.end(),
	       baseline, bsBegin, bsEnd);

      Rectangle maxRect;
      maxRect.setNumAnomalies(0, 0);
      maxRect.setNumPoints(0, 0);
      int aTotal = asEnd - asBegin;
      int bTotal = bsEnd - bsBegin;
      double eps = 1.0 / static_cast<double>(end - begin);
      double maxScan = 0;
      auto earlyTerm = [&anomalies, &baseline,
			eps, rho,
			aTotal, bTotal](int i, int j,
					int k, int m,
					double ms) {
      	int aCount = computeCell(anomalies, i, j, k, m);
	int bCount = computeCell(baseline, i, j, k, m);
	double mr = (double) aCount / aTotal;
	double br = (double) bCount / bTotal;
	return ms < maxSub(mr, br, eps, rho);
      };

      for (int i = 0; i < g - 1 && earlyTerm(i, 0, g - 1, g - 1, maxScan); i++) { // left
	for (int j = 0; j < g - 1 && earlyTerm(i, j, g - 1, g - 1, maxScan); j++) { // top
	  for (int k = g - 1; k > i && earlyTerm(i, j, k, g - 1, maxScan); k--) { // right
	    for (int m = g - 1; m > j && earlyTerm(i, j, k, m, maxScan); m--) { // bottom
	      int aCount = computeCell(anomalies, i, j, k, m);
	      int bCount = computeCell(baseline, i, j, k, m);
	      double mr = (double) aCount / aTotal;
	      double br = (double) bCount / bTotal;
	      if (rhoInRange(mr, br, rho, eps)) {
		Rectangle rect(x_pos[i], x_pos[k], y_pos[j], y_pos[m]);
		rect.setNumAnomalies(aCount, aTotal);
		rect.setNumPoints(bCount, bTotal);
		double newSt = rect.statistic();
		if (newSt >= maxScan) {
		  maxRect = rect;
		  maxScan = newSt;
		}
	      }
	    }
	  }
	}
      }
      return maxRect;
    }
  };
}

#endif
