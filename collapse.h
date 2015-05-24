#ifndef COLLAPSE_H
#define COLLAPSE_H

#include <iostream>

#include "metric_space.h"
#include "simplicial_set.h"

using namespace std;

/* A way of collapsing a simplicial cx, as well as
 * collapsing a cx in a filtration in a way compatible
 * with the rest of the collapsed cxs in the filtration
 */
class CollapseStrategy {
public:
  virtual SimplicialSetPtr collapse_fully(CliqueGraphPtr cgp) = 0;
  virtual SimplicialSetPtr collapse_partial(CliqueGraphPtr cgp,
                                            SimplicialSetPtr previous) = 0;
};

/* method to construct a collapsed filtered cx given appropriate
 * data, a clique factory, and collapsing strategy.
 * epsilons is sorted from smallest to largest
 */
template <typename CliqueFactory, typename CollapseStrategy>
FilteredSimplicialSetPtr
collapsed_filtered_cx(FiniteMetricSpace& fms,
                      vector<coord_t>& epsilons,
                      int k,
                      CliqueFactory clique_factory,
                      CollapseStrategy& collapser) {
  auto fss = make_shared<FilteredSimplicialSet>(k);
  auto max_epsilon = *epsilons.rbegin();
  WeightedNeighborhoodGraph ng(fms, max_epsilon);

  for (auto e = epsilons.rbegin(); e != epsilons.rend(); ++e) {
    auto epsilon = *e;
    auto cgp = ng.get_clique_graph(clique_factory, epsilon);
    if (e == epsilons.rbegin()) {
      fss->add_simplicial_set(epsilon, collapser.collapse_fully(cgp));
    } else {
      auto ss_previous = fss->last_simplical_set();
      auto cur = collapser.collapse_partial(cgp, ss_previous);
      fss->add_simplicial_set(epsilon, cur);
    }
  }
  return fss;
}

/* collapse NO simplices
 */
class TrivialCollapseStrategy : public CollapseStrategy {
public:
  SimplicialSetPtr collapse_fully(CliqueGraphPtr cgp) {
    vector<Simplex> empty;
    vector<Simplex> vertices;
    vertices.reserve(cgp->n_vertices());
    cgp->get_vertices(back_inserter(vertices));
    return make_shared<SimplicialSet>(vertices, empty);
  }

  SimplicialSetPtr collapse_partial(CliqueGraphPtr cgp,
                                    SimplicialSetPtr previous) {
    return collapse_fully(cgp);
  }
};

/* collapse sxs by choosing appropriate nodes in the clique
 * graph.  Do no homology calculations.
 */
class NaiveCollapseStrategy : public CollapseStrategy {
protected:
  template <typename Collection>
  bool zero_or_one_in(Collection& c1, Collection& c2) {
    int seen = 0;
    for (auto i = c1.begin(); i != c1.end(); ++i) {
      if (c2.find(*i) != c2.end())
        seen++;
      if (seen > 1)
        return false;
    }
    return true;
  }

  template <typename OutputIterator>
  void get_collapsable_vertices(CliqueGraphPtr cgp,
                                OutputIterator out,
                                SimplicialSetPtr ss) {
    set<Simplex> collapsed;
    vector<Simplex> verts;
    set<Simplex> nbrs;

    verts.reserve(cgp->n_vertices());

    // New way, doesn't seem to work as well...
    //cgp->vertices_by_dimension(back_inserter(verts));

    // Old way; seems to work better in terms
    // of creating sx sets w/ fewer sxs...
    Simplex max_vert = cgp->maximal_vertex();
    cgp->get_vertices(inserter(verts, verts.begin()));

    if (max_vert.dimension() == 0)
      return; // no collapsing can be performed

    verts.erase(remove(verts.begin(), verts.end(), max_vert),
                verts.end());
    *(out++) = max_vert;
    collapsed.insert(max_vert);
    // end old way

    for (auto i = verts.begin(); i != verts.end(); ++i) {
      if (i->dimension() == 0) // don't collapse points
        continue;
      nbrs.clear();
      cgp->get_neighbors(*i, inserter(nbrs, nbrs.begin()));
      if (zero_or_one_in(nbrs, collapsed)) {
        if (ss && !ss->is_collapsed(*i))
          continue;
        collapsed.insert(*i);
        *(out++) = *i;
      }
    }
  }

  SimplicialSetPtr collapse(CliqueGraphPtr cgp, SimplicialSetPtr ss) {
    vector<Simplex> all_vertices;
    vector<Simplex> collapsed_vertices;
    vector<Simplex> maximal_vertices;

    all_vertices.reserve(cgp->n_vertices());
    cgp->get_vertices(back_inserter(all_vertices));

    get_collapsable_vertices(cgp, back_inserter(collapsed_vertices), ss);

    sort(all_vertices.begin(), all_vertices.end());
    sort(collapsed_vertices.begin(), collapsed_vertices.end());
    set_difference(all_vertices.begin(), all_vertices.end(),
                   collapsed_vertices.begin(), collapsed_vertices.end(),
                   back_inserter(maximal_vertices));

    return make_shared<SimplicialSet>(maximal_vertices,
                                      collapsed_vertices);
  }
public:
  SimplicialSetPtr collapse_fully(CliqueGraphPtr cgp) {
    SimplicialSetPtr ss;
    return collapse(cgp, ss);
  }

  SimplicialSetPtr collapse_partial(CliqueGraphPtr cgp,
                                    SimplicialSetPtr previous) {
    return collapse(cgp, previous);
  }

};

/****************************************/

template <typename Vertex, typename Iterator>
class DfsCollapseVisitor : public default_dfs_visitor {
protected:
  SimplicialSetPtr ss;
  Iterator out;

  std::unordered_map<Vertex, int> collapsed_nbrs_map;
  Vertex last_vertex;
public:
  DfsCollapseVisitor(SimplicialSetPtr ss, Iterator out) : ss(ss), out(out) {
  }

  template <typename Graph>
  void discover_vertex(Vertex u, const Graph & g) {
    if (collapsed_nbrs_map[u] < 2) {
      Simplex s = g[u].simplex;
      if (ss && !ss->is_collapsed(s))
        return;
      if (s.dimension() > 0) {
        *out++ = s;
        typedef typename graph_traits<Graph>::adjacency_iterator
          AdjacencyIterator;
        AdjacencyIterator vi, vi_end;
        for (tie(vi, vi_end) = adjacent_vertices(u, g);
             vi != vi_end; ++vi) {
          collapsed_nbrs_map[*vi] += 1;
        }
      }
    }
  }
};

class HomologyCollapseVisitorStrategy : public CollapseStrategy {
protected:
  template <typename OutputIterator>
  void get_collapsable_vertices(CliqueGraphPtr cgp,
                                OutputIterator out,
                                SimplicialSetPtr ss) {
    DfsCollapseVisitor<CliqueGraph::vertex_t, OutputIterator> v(ss, out);
    cgp->dfs(v);
  }

  SimplicialSetPtr collapse(CliqueGraphPtr cgp, SimplicialSetPtr ss) {
    vector<Simplex> all_vertices;
    vector<Simplex> collapsed_vertices;
    vector<Simplex> maximal_vertices;

    all_vertices.reserve(cgp->n_vertices());
    cgp->get_vertices(back_inserter(all_vertices));

    get_collapsable_vertices(cgp, back_inserter(collapsed_vertices), ss);

    sort(all_vertices.begin(), all_vertices.end());
    sort(collapsed_vertices.begin(), collapsed_vertices.end());
    set_difference(all_vertices.begin(), all_vertices.end(),
                   collapsed_vertices.begin(), collapsed_vertices.end(),
                   back_inserter(maximal_vertices));

    // XXX this is horrible and needs to be rewritten

    DegenerateSimplices ds(collapsed_vertices.begin(),
                           collapsed_vertices.end());

    for (auto mv = maximal_vertices.begin();
         mv != maximal_vertices.end(); ++mv) {
      if (ds.is_acyclic(*mv)) {
        ds.add_simplex(*mv);
        collapsed_vertices.push_back(*mv);
      }
    }

    maximal_vertices.clear();
    sort(collapsed_vertices.begin(), collapsed_vertices.end());
    set_difference(all_vertices.begin(), all_vertices.end(),
                   collapsed_vertices.begin(), collapsed_vertices.end(),
                   back_inserter(maximal_vertices));

    return make_shared<SimplicialSet>(maximal_vertices,
                                      collapsed_vertices);

  }
public:
  SimplicialSetPtr collapse_fully(CliqueGraphPtr cgp) {
    SimplicialSetPtr ss;
    return collapse(cgp, ss);
  }

  SimplicialSetPtr collapse_partial(CliqueGraphPtr cgp,
                                    SimplicialSetPtr previous) {
    return collapse(cgp, previous);
  }
};

/**************************************/

/*  collapse using combinatorial strategy, then
 *  look at the top percent of maximal simplices by dimension,
 *  try to collapse those
 */
class PartialCollapseStrategy : public CollapseStrategy {
public:
  PartialCollapseStrategy(float percent)
    : percent(percent) {}
protected:
  double percent;

  template <typename OutputIterator>
  void get_collapsable_vertices(CliqueGraphPtr cgp,
                                OutputIterator out,
                                SimplicialSetPtr ss) {
    DfsCollapseVisitor<CliqueGraph::vertex_t, OutputIterator> v(ss, out);
    cgp->dfs(v);
  }

  SimplicialSetPtr collapse(CliqueGraphPtr cgp, SimplicialSetPtr ss) {
    vector<Simplex> all_vertices;
    vector<Simplex> collapsed_vertices;
    vector<Simplex> maximal_vertices;

    all_vertices.reserve(cgp->n_vertices());
    cgp->get_vertices(back_inserter(all_vertices));

    get_collapsable_vertices(cgp, back_inserter(collapsed_vertices), ss);

    sort(all_vertices.begin(), all_vertices.end());
    sort(collapsed_vertices.begin(), collapsed_vertices.end());
    set_difference(all_vertices.begin(), all_vertices.end(),
                   collapsed_vertices.begin(), collapsed_vertices.end(),
                   back_inserter(maximal_vertices));

    vector<Simplex> top_dimension;
    vector<Simplex> try_to_collapse;
    int n = cgp->n_vertices() * percent;
    cout << " trying to collapse " << n << " / " <<
      cgp->n_vertices() << endl;

    cgp->maximal_vertices(n, back_inserter(top_dimension));
    sort(top_dimension.begin(), top_dimension.end());
    cout <<  " top_dim size " << top_dimension.size() << endl;
    set_difference(top_dimension.begin(), top_dimension.end(),
                   collapsed_vertices.begin(),
                   collapsed_vertices.end(),
                   back_inserter(try_to_collapse));

    cout << " actually trying " << try_to_collapse.size() << endl;

    DegenerateSimplices ds(collapsed_vertices.begin(),
                           collapsed_vertices.end());

    for (auto mv = try_to_collapse.begin();
         mv != try_to_collapse.end(); ++mv) {
      if (ds.is_acyclic(*mv)) {
        ds.add_simplex(*mv);
        collapsed_vertices.push_back(*mv);
      }
    }

    maximal_vertices.clear();
    sort(collapsed_vertices.begin(), collapsed_vertices.end());
    set_difference(all_vertices.begin(), all_vertices.end(),
                   collapsed_vertices.begin(), collapsed_vertices.end(),
                   back_inserter(maximal_vertices));

    return make_shared<SimplicialSet>(maximal_vertices,
                                      collapsed_vertices);

  }
public:
  SimplicialSetPtr collapse_fully(CliqueGraphPtr cgp) {
    SimplicialSetPtr ss;
    return collapse(cgp, ss);
  }

  SimplicialSetPtr collapse_partial(CliqueGraphPtr cgp,
                                    SimplicialSetPtr previous) {
    return collapse(cgp, previous);
  }
};

/* This is old; walks the clique graph explicitly.  We should
 * probably do it using the visitor as above.
 */
class HomologyCollapseStrategy : public CollapseStrategy {
protected:
  template <typename Collection>
  bool zero_or_one_in(Collection& c1, Collection& c2) {
    int seen = 0;
    for (auto i = c1.begin(); i != c1.end(); ++i) {
      if (c2.find(*i) != c2.end())
        seen++;
      if (seen > 1)
        return false;
    }
    return true;
  }

  // XXXX ARG
  template <typename OutputIterator>
  void get_collapsable_vertices(CliqueGraphPtr cgp,
                                OutputIterator out,
                                SimplicialSetPtr ss) {
    set<Simplex> visited;
    set<Simplex> collapsed;
    set<Simplex> neighbors;
    set<Simplex> all_vertices;

    stack<Simplex> to_check;

    cgp->get_vertices(inserter(all_vertices, all_vertices.begin()));
    Simplex max_vertex = cgp->maximal_vertex();
    to_check.push(max_vertex);

    while (!to_check.empty()) {
      Simplex s = to_check.top();
      to_check.pop();
      visited.insert(s);

      if (s.dimension() == 0) // don't collapse points
        continue;

      cgp->get_neighbors(s, inserter(neighbors, neighbors.begin()));
      if (zero_or_one_in(neighbors, collapsed)) {
        if (ss && !ss->is_collapsed(s))
          continue;
        *(out++) = s;
        collapsed.insert(s);
      }
      for (auto n = neighbors.begin(); n != neighbors.end(); ++n)
        if (visited.count(*n) == 0)
          to_check.push(*n);
      if (to_check.empty()) { // try to add vertex from another connected cpt
        // XXX probably a better way to do this
        if (visited.size() == all_vertices.size())
          return;
        for (auto v = all_vertices.begin(); v != all_vertices.end(); ++v) {
          if (visited.count(*v) == 0) {
            to_check.push(*v);
            break;
          }
        }
      }
    }
  }


  SimplicialSetPtr collapse(CliqueGraphPtr cgp, SimplicialSetPtr ss) {
    vector<Simplex> all_vertices;
    vector<Simplex> collapsed_vertices;
    vector<Simplex> maximal_vertices;

    all_vertices.reserve(cgp->n_vertices());
    cgp->get_vertices(back_inserter(all_vertices));

    get_collapsable_vertices(cgp, back_inserter(collapsed_vertices), ss);

    sort(all_vertices.begin(), all_vertices.end());
    sort(collapsed_vertices.begin(), collapsed_vertices.end());
    set_difference(all_vertices.begin(), all_vertices.end(),
                   collapsed_vertices.begin(), collapsed_vertices.end(),
                   back_inserter(maximal_vertices));

    return make_shared<SimplicialSet>(maximal_vertices,
                                      collapsed_vertices);

  }
public:
  SimplicialSetPtr collapse_fully(CliqueGraphPtr cgp) {
    SimplicialSetPtr ss;
    return collapse(cgp, ss);
  }

  SimplicialSetPtr collapse_partial(CliqueGraphPtr cgp,
                                    SimplicialSetPtr previous) {
    return collapse(cgp, previous);
  }

};

#endif
