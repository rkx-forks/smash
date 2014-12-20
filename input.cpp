#include <sstream>
#include <iostream>
#include <iterator>
#include <cassert>
#include <boost/algorithm/string.hpp>

#include "metric_space.h"

using namespace std;

void read_points(std::istream &str, point_list &pl) {
  string line;
  int dimension = -1;
  
  while (getline(str, line)) {
    stringstream ss(line);
    point p((istream_iterator<coord_t>(ss)), istream_iterator<coord_t>());
    if (dimension == -1) {
      dimension = p.size();
    } else if (dimension != (int)p.size()) {
      cout << "expected " << dimension << ", got " << p.size() << endl;
      assert(false);
    }
    pl.push_back(p);
    //    boost::split(nums, line, boost::is_any_of("\t "));
  }
}
