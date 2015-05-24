#ifndef CLIQUE_H
#define CLIQUE_H

#include <vector>
#include <set>
#include <map>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/lookup_edge.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/connected_components.hpp>

#include "simplex.h"

using namespace boost;

// annoying: can't really get to vertex type w/o template parameter...
/*typedef long unsigned int boost_vertex_t;

typedef std::vector<boost_vertex_t> clique_graph_vertex_t;
typedef std::vector<clique_graph_vertex_t> clique_vertex_container_t;*/


/*
  a graph whose vertices are simplices and edges represent
  intersection between simplices (Maybe should be called "Dual Graph")
 */
struct CGVertexInfo {
  Simplex simplex;
};

class CliqueGraph {
public:
  typedef adjacency_list<vecS, vecS, undirectedS,
                         CGVertexInfo> clique_graph_t;
  typedef graph_traits<clique_graph_t>::vertex_descriptor vertex_t;


protected:
  typedef std::unordered_map<Simplex, vertex_t, SimplexHasher> vertex_map_t;

  clique_graph_t g;
  vertex_map_t vertex_map;
 public:
  template <typename Iterator>
  CliqueGraph(Iterator begin, Iterator end) {
    for (auto i = begin; i != end; ++i) {
      vertex_t v = add_vertex(g);
      vertex_map[*i] = v;
      g[v].simplex = *i;
    }

    for (auto i = vertex_map.begin(); i != vertex_map.end(); ++i) {
      auto j = i;
      j++;
      Simplex s1 = i->first;
      for (; j != vertex_map.end(); ++j) {
        Simplex s2 = j->first;
        if (s1.has_intersection(s2))
          add_edge(i->second, j->second, g);
      }
    }

    assert(!is_degenerate());
  }

  bool is_degenerate() { // check to make sure no empty simplices
    for (auto vi = vertices(g); vi.first != vi.second; ++vi.first) {
      Simplex s = g[*vi.first].simplex;
      if (s.dimension() < 0)
        return true;
    }
    return false;
  }

  bool has_vertex(Simplex& s) {
    return vertex_map.count(s) > 0;
  }

  bool has_edge(Simplex& s1, Simplex& s2) {
    bool found;
    graph_traits<clique_graph_t>::edge_descriptor e;
    tie(e, found) = edge(vertex_map[s1], vertex_map[s2], g);
    return found;
  }

  vertex_map_t::size_type n_vertices(void) {
    return vertex_map.size();
  }

  template <typename OutputIterator>
  void get_vertices(OutputIterator out) {
    for (auto i = vertex_map.begin(); i != vertex_map.end(); ++i)
      *(out++) = i->first;
  }

  template <typename OutputIterator>
  void get_neighbors(Simplex& s, OutputIterator out) {
    assert(vertex_map.find(s) != vertex_map.end());
    clique_graph_t::adjacency_iterator i, i_end;
    tie(i, i_end) = adjacent_vertices(vertex_map[s], g);
    for (; i != i_end; ++i)
      *(out++) = g[*i].simplex;
  }

  template <typename OutputIterator>
  void vertices_by_dimension(OutputIterator out) {
    vector<Simplex> sxs;
    for (auto i = vertex_map.begin(); i != vertex_map.end(); ++i)
      sxs.push_back(i->first);
    sort(sxs.begin(), sxs.end(), [](const Simplex& s1, const Simplex& s2) {
        return s1.dimension() > s2.dimension();
      });
    cout << "Sxs by dim: " << endl;
    for (int i = 0; i < 5; ++i) {
      cout << sxs[i] << endl;
    }
    for (const auto& sx : sxs) {
      *(out++) = vertex_map[sx];
    }
  }

  Simplex maximal_vertex(void) {
    int max_dim = -1;
    Simplex max_sx;
    for (auto i = vertex_map.begin(); i != vertex_map.end(); ++i) {
      if (i->first.dimension() > max_dim) {
        max_sx = i->first;
        max_dim = max_sx.dimension();
      }
    }
    return max_sx;
  }

  /* return  n simplexes sorted by
     decreasing dimension
   */
  template <typename Iterator>
  void maximal_vertices(size_t n, Iterator out) {
    vector<Simplex> sxs;
    transform(vertex_map.begin(), vertex_map.end(),
              back_inserter(sxs),
              [&](pair<const Simplex, vertex_t>& p) {
                return p.first;
              });
    sort(sxs.begin(), sxs.end(),
         [&](const Simplex& s1, const Simplex& s2) {
           return s1.dimension() > s2.dimension();
         });
    assert(n < sxs.size());
    copy(sxs.begin(), sxs.begin() + n, out);
  }

  void all_simplices(int max_dimension, set<Simplex>& output) {
    for (auto i = vertex_map.begin(); i != vertex_map.end(); ++i) {
      Simplex s = (*i).first;
      s.all_subsimplices(max_dimension, inserter(output, output.end()));
    }
  }

  template <typename Visitor>
  void dfs(Visitor v) {
    depth_first_search(g, visitor(v));
  }

  void print_verts() {
    std::cout << "CliqueGraph verts ";
    for (auto i = vertex_map.begin(); i != vertex_map.end(); ++i) {
      std::cout << i->first << endl;
    }
  }

  void print() {
    std::cout << "CliqueGraph: " << endl;
    std::cout << "vertices:" << endl;
    for (auto i = vertex_map.begin(); i != vertex_map.end(); ++i)
      std::cout << i->first << endl;

    clique_graph_t::edge_iterator e, e_end;
    for (tie(e, e_end) = edges(g); e != e_end; ++e) {
      std::cout << g[source(*e, g)].simplex << " - "
                << g[target(*e, g)].simplex << endl;
    }
  }

  /*  template <typename Iterator>
  void connected_component_vertices(Iterator out) {
    vector<int> component(num_vertices(g));
    vector<int>::size_type i;
    int num = connected_components(g, &component[0]);
    vector< set<Simplex::vertex_t> > verts(num);

    for (i = 0; i < component.size(); ++i) {
      copy(g[i].simplex.begin(), g[i].simplex.end(),
           inserter(verts[component[i]], verts[component[i]].begin()));
    }
    for (int c = 0; c < num; ++c) {
      *out++ = Simplex(verts[c].begin(), verts[c].end());
    }
    }*/
};

typedef std::shared_ptr<CliqueGraph> CliqueGraphPtr;


#endif
