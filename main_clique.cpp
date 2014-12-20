#include <iostream>
#include <boost/lexical_cast.hpp>

#include "metric_space.h"
#include "networkx_clique.h"
#include "incremental_vr.h"
#include "input.h"

int main(int argc, char **argv) {
  SimplicialComplex cx;
  point_list pts;

  if (argc != 3) {
    std::cout << "use: " << argv[0] << " epsilon k < [input data]" << std::endl;
    return 1;
  }

  coord_t epsilon = boost::lexical_cast<coord_t>(argv[1]);
  int k = boost::lexical_cast<int>(argv[2]);

  read_points(std::cin, pts);
  FiniteMetricSpace fms(2, pts);

  NeighborhoodGraphPtr ngp = fms.get_neighborhood_graph(epsilon*epsilon);
  
  NetworkxCliqueFactory clique_factory;
  //BKCliqueFactory clique_factory;
  CliqueGraphPtr cgp = ngp->get_clique_graph(clique_factory);

  //ngp->print();
  //  cgp->print();

  //  ngp->print();
  //std::cout << "clique vertices: " << std::endl;
  //  cgp->print_verts();
  //  std::cout << "end clique verts " << std::endl;
  set<Simplex> sxs;
  cgp->all_simplices(k, sxs);
   
 for (set<Simplex>::iterator i = sxs.begin();
       i != sxs.end(); ++i) {
    std::cout << *i << std::endl;
  }
}
