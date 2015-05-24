#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include "metric_space.h"
#include "input.h"
#include "networkx_clique.h"
#include "collapse.h"
#include "simplicial_set.h"

namespace po = boost::program_options;

int main(int argc, char **argv) {
  bool verbose;
  int k;
  vector<coord_t> epsilons;
  string collapse_strategy;

  po::options_description opts("Allowed options");
  opts.add_options()
    ("help", "produce help message")
    ("verbose,v", po::value<bool>(&verbose)->zero_tokens())
    ("strategy,s", po::value<string>(&collapse_strategy)
        ->default_value("naive"))
    ;

  po::options_description hidden_opts("Hidden opts");
  hidden_opts.add_options()
    ("k", po::value<int>(&k), "(positional args, don't use)")
    ("epsilons", po::value< vector<coord_t> >(&epsilons))
    ;

  po::options_description all_options("All Options");
  all_options.add(opts).add(hidden_opts);

  po::variables_map vm;
  po::positional_options_description p;
  p.add("k", 1);
  p.add("epsilons", -1);
  auto parsed = po::command_line_parser(argc, argv)
    .options(all_options).positional(p).run();
  po::store(parsed, vm);
  po::notify(vm);

  if (vm.count("help") || !epsilons.size()) {
    cout << "Usage: " << argv[0] << " k e1 e2 ..." << endl;
    cout << opts << "\n";
    return 1;
  }

  unique_ptr<CollapseStrategy> collapse_strategy_p;
  if (collapse_strategy == "naive" || collapse_strategy == "n") {
    collapse_strategy_p.reset(new NaiveCollapseStrategy);
    cout << "Using naive collapse strategy" << endl;
  } else if (collapse_strategy == "partial" || collapse_strategy == "p") {
    collapse_strategy_p.reset(new PartialCollapseStrategy(0.25));
    cout << "Using partial collapse strategy" << endl;
  } else if (collapse_strategy == "homology" || collapse_strategy == "h") {
    collapse_strategy_p.reset(new HomologyCollapseVisitorStrategy);
    cout << "Using homology collapse strategy" << endl;
  } else if (collapse_strategy == "trivial" || collapse_strategy == "t") {
    collapse_strategy_p.reset(new TrivialCollapseStrategy);
    cout << "Using trivial collapse strategy" << endl;
  } else {
    cout << "invalid collapse strategy '" << collapse_strategy << "'" << endl;
    cout << "valid: naive | partial | homology | trivial" << endl;
    return 1;
  }

  point_list pts;
  read_points(std::cin, pts);
  FiniteMetricSpace fms(pts);

  FilteredSimplicialSetPtr fss =
    collapsed_filtered_cx(fms, epsilons, k, NetworkxCliqueFactory(),
                          *collapse_strategy_p);

  if (verbose)  {
    fss->compute_p_intervals();
    fss->print_info();
    //fss->print_induced_maps();
    //fss->print_p_intervals();
    cout << "(Imagine p-intervals printed here)" << endl;
  }
}
