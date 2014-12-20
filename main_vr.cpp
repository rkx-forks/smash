#include <iostream>
#include <boost/lexical_cast.hpp>

#include "metric_space.h"
#include "incremental_vr2.h"
#include "input.h"

int main(int argc, char **argv) {
  //SimplicialComplex cx;
  point_list pts;

  if (argc != 3) {
    std::cout << "use: " << argv[0] << " epsilon k < [input data]" << std::endl;
    return 1;
  }

  coord_t epsilon = boost::lexical_cast<coord_t>(argv[1]);
  int k = boost::lexical_cast<int>(argv[2]);

  read_points(std::cin, pts);
  FiniteMetricSpace fms(pts);
  WeightedNeighborhoodGraph ng(fms, epsilon);

  //cout << ng << endl;
  TestingVrVisitor visitor;
  incremental_vr(ng, k, visitor);
  
  visitor.print();

  //cx.incremental_vr(fms, epsilon, k);
  //cx.print();
}
