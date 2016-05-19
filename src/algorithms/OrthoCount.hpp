#ifndef ORTHOCOUNT
#define ORTHOCOUNT
#include <vector>
#include <algorithm>
#include <tuple>
#include <deque>

#include "framework/Point.hpp"

namespace anomaly {
  
  using namespace std;

  /*
  struct {
    bool operator() (Point const& pt1, Point const& pt2) {
      return pt1.getY() != pt2.getY() ? pt1.getY() < pt2.getY() : pt1.getX() < pt2.getX();
    }
  } cmp_yx;
    
  struct {
    bool operator() (Point const& pt1, Point const& pt2) {
      return pt1.getX() != pt2.getX() ? pt1.getX() < pt2.getX() : pt1.getY() < pt2.getY();
    }
  } cmp_xy;
  */
  struct {
    bool operator() (Point const& pt1, Point const& pt2) {
      return pt1.getY() < pt2.getY();
    }
  } cmp_yx;
    
  struct {
    bool operator() (Point const& pt1, Point const& pt2) {
      return pt1.getX() < pt2.getX();
    }
  } cmp_xy;

  class OrthoCount {
    using pt_t = Point;

    struct RangeNode {
      pt_t median_pt;
      shared_ptr<RangeNode> right;
      shared_ptr<RangeNode> left;
      vector<long> left_map;
      vector<long> right_map;
      vector<Point> pt_set;
	
      long getRightIndex(long const& index) {
	return right_map[index];
      }

      long getLeftIndex(long const& index) {
	return left_map[index];
      }

      auto getPoint(long const& index) -> decltype(pt_set.begin()) {
	return pt_set.begin() + index;
      }
    };
    
    shared_ptr<RangeNode> root;
    vector<pt_t> sorted_xy;
    
  public:

    OrthoCount() {}
    
    template<typename T>
    OrthoCount(T begin, T end) :
      sorted_xy(begin, end) {
      /*
	This will modify its template arguments.
      */
      //Sort these points for the recursive descent into the tree.
      sort(sorted_xy.begin(), sorted_xy.end(), cmp_xy);
      vector<tuple<T, T, vector<pt_t>, shared_ptr<RangeNode>*>> nodeStack;
      nodeStack.emplace_back(begin, end, sorted_xy, &root);
      
      while (nodeStack.size() > 0) {
	auto curr_node = nodeStack.back();
	nodeStack.pop_back();
	auto& curr_begin_yx = get<0>(curr_node);
	auto& curr_end_yx = get<1>(curr_node);
	auto& curr_pts_xy = get<2>(curr_node);
	auto curr_address = get<3>(curr_node);
	if (curr_end_yx != curr_begin_yx) {
	  *curr_address = make_shared<RangeNode>();
	  (*curr_address)->pt_set = curr_pts_xy;

	  if (curr_end_yx - curr_begin_yx > 1) {
	    //Split the points based on the median.
	    nth_element(curr_begin_yx, curr_begin_yx + (curr_end_yx - curr_begin_yx) / 2, curr_end_yx, cmp_yx);

	    //Create tuples for the stack
	    nodeStack.emplace_back(curr_begin_yx, curr_begin_yx + (curr_end_yx - curr_begin_yx) / 2,
				   vector<pt_t>(curr_begin_yx, curr_begin_yx + (curr_end_yx - curr_begin_yx) / 2),
				   &((*curr_address)->left));
	    auto lpts_begin = get<2>(nodeStack.back()).begin();
	    auto lpts_end = get<2>(nodeStack.back()).end();
	    nodeStack.emplace_back(curr_begin_yx + (curr_end_yx - curr_begin_yx) / 2, curr_end_yx,
				   vector<pt_t>(curr_begin_yx + (curr_end_yx - curr_begin_yx) / 2 + 1, curr_end_yx),
				   &((*curr_address)->right));
	    auto rpts_begin = get<2>(nodeStack.back()).begin();
	    auto rpts_end = get<2>(nodeStack.back()).end();
	    //Now sort the left and right vectors and make a mapping
	    sort(lpts_begin, lpts_end, cmp_xy); 
	    sort(rpts_begin, rpts_end, cmp_xy);
	    (*curr_address)->left_map.resize(end - begin);
	    (*curr_address)->right_map.resize(end - begin);
	    //Now generate the parent child mapping.
	    long r = 0;
	    long l = 0;
	    for (long i = 0; i < curr_pts_xy.size(); i++) {
	      (*curr_address)->right_map[i] = r;
	      (*curr_address)->left_map[i] = l;
	      if (r < rpts_end - rpts_begin && rpts_begin[r] == curr_pts_xy[i]) {
		r++;
	      } else {
		l++;
	      }
	    }
	  }
	  (*curr_address)->median_pt = *(curr_begin_yx + (curr_end_yx - curr_begin_yx) / 2);
	}	
      }
    }

    template <typename F>
    void foldRange(pt_t const& upper, pt_t const& lower, F f) {
      auto& left_node = root;
      auto& right_node = root;
      //Run upper and lower bound to get start indices into root
      // Iterate down through the children following the maps.
      // check median then follow index

      auto it_up = upper_bound(sorted_xy.begin(), sorted_xy.end(), upper, cmp_xy);
      auto it_lower = upper_bound(sorted_xy.begin(), sorted_xy.end(), lower, cmp_xy);
      long left_lower_index = it_lower - sorted_xy.begin(); 
      long right_upper_index = it_up - sorted_xy.begin();
      long left_upper_index = it_up - sorted_xy.begin(); 
      long right_lower_index = it_lower - sorted_xy.begin();

      while (left_node || right_node) {
	if (left_node) {
	  if (cmp_xy(lower, left_node->median_pt)) {
	    if (left_node != right_node) {
	      f(left_node->right,
		left_node->getRightIndex(left_lower_index),
		left_node->getRightIndex(left_upper_index));
	    }
	    left_upper_index = left_node->getLeftIndex(left_upper_index);
	    left_lower_index = left_node->getLeftIndex(left_lower_index);
	    left_node = left_node->left;
	  } else {
	    left_upper_index = left_node->getRightIndex(left_upper_index);
	    left_lower_index = left_node->getRightIndex(left_lower_index);
	    left_node = left_node->right;
	  }
	}
	if (right_node) {
	  if (cmp_xy(right_node->median_pt, upper)) {
	    if (left_node != right_node) {
	      f(right_node->left,
		right_node->getLeftIndex(right_lower_index),
		right_node->getLeftIndex(right_upper_index));
	    }
	  } else {
	    right_upper_index = right_node->getLeftIndex(right_upper_index);
	    right_lower_index = right_node->getLeftIndex(right_lower_index);
	    right_node = right_node->left;
	  }
	}
      }
    }

    long countRange(pt_t const& upper, pt_t const& lower) {
      long count = 0;
      foldRange(upper, lower, [&count](shared_ptr<RangeNode> const& node, long lower, long upper){
	  count += upper - lower;
	});
      return count;
    }

    template<typename T>
    void accumRange(T begin, T end, pt_t const& upper, pt_t const& lower) {
      foldRange(upper, lower, [&begin, &end](shared_ptr<RangeNode> const& node, long lower, long upper){
	  if (node) {
	    auto ed = node->getPoint(upper);
	    auto bg = node->getPoint(lower);
	    while (bg != ed) {
	      *begin = *bg;
	      begin++;
	      bg++;
	    }
	  }
	});
    }

    friend ostream& operator<<(ostream& os, OrthoCount& oc);
  };
  
  ostream& operator<<(ostream& os, OrthoCount& oc) {
    deque<tuple<shared_ptr<OrthoCount::RangeNode>, int>> nodeQueue;
    nodeQueue.push_back(make_tuple(oc.root, 0));
    int currDepth = 0;
    while (nodeQueue.size() > 0) {
      auto first_el = nodeQueue.front();
      if (currDepth < get<1>(first_el)) {
	currDepth += 1;
	os << endl;
      }

      os << "{";
      for (auto& pt : get<0>(first_el)->pt_set) {
	os << pt << " , ";
      }
      os << "} | " ;
      if (get<0>(first_el)->right) {
	nodeQueue.push_back(make_tuple(get<0>(first_el)->right, get<1>(first_el) + 1));
      }
      if (get<0>(first_el)->left) {
	nodeQueue.push_back(make_tuple(get<0>(first_el)->left, get<1>(first_el) + 1));
      }
      
      nodeQueue.pop_front();
	
    }
    return os;
  }

}
#endif
