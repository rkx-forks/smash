#include "metric_space.h"
#include "clique.h"

struct CliqueVisitor {
  template <typename Clique, typename Graph>
  void clique(Clique& p, Graph& g) {
    Simplex s(p.begin(), p.end());
    simplices.push_back(s);
  }

  vector<Simplex> simplices;
};

template <typename Graph, typename Visitor>
void bron_kerbosch_all_cliques(const Graph& g, Visitor& visitor);

struct BKCliqueFactory {
  template <typename Graph>
  CliqueGraphPtr make_clique_graph(Graph& graph) {
    CliqueVisitor visitor;
    bron_kerbosch_all_cliques(graph, visitor);
    return make_shared<CliqueGraph>(visitor.simplices.begin(),
                                    visitor.simplices.end());
  }
};


/* old BK stuff:
 */
template <typename Graph, typename Container>
  void filter_unconnected_vertices(const Graph& g, 
                                   typename graph_traits<Graph>::vertex_descriptor v, 
                                   const Container& in,
                                   Container& out) {
  typename Container::const_iterator i, end = in.end();
  for (i = in.begin(); i != end; ++i) {
    if (lookup_edge(v, *i, g).second) {
      out.push_back(*i);
    }
  }
}

template <typename Graph, typename Clique, typename Container, 
  typename Visitor>
  void extend_clique(const Graph& g, Clique& clique, Container& cands,
                     Container& nots, Visitor& visitor,
                     std::size_t min) {
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  
  // vertex in nots connected to all vertices in cands?
  typename Container::iterator ni, nend = nots.end();
  typename Container::iterator ci, cend = cands.end();
  for (ni = nots.begin(); ni != nend; ++ni) {
    for (ci = cands.begin(); ci != cend; ++ci) {
      if (!lookup_edge(*ni, *ci, g).second) break;
    }
    if (ci == cend) break;
  }
  if (ni != nend) return;

  // loop through cands
  typename Container::iterator i, j;
  for (i = cands.begin(); i != cands.end();) {
    Vertex candidate = *i;
    clique.push_back(candidate);
    i = cands.erase(i);
    Container new_cands, new_nots;
    filter_unconnected_vertices(g, candidate, cands, new_cands);
    filter_unconnected_vertices(g, candidate, nots, new_nots);
    
    if (new_cands.empty() && new_nots.empty()) {
      if (clique.size() >= min) 
        visitor.clique(clique, g);
    } else {
      extend_clique(g, clique, new_cands, new_nots, visitor, min);
    }
    nots.push_back(candidate);
    clique.pop_back();
  }
}


template <typename Graph, typename Visitor>
void bron_kerbosch_all_cliques(const Graph& g, 
                               Visitor& visitor) {
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertex_iterator VertexIterator;
  typedef std::vector<Vertex> VertexSet;
  typedef std::vector<Vertex> Clique;
 
  VertexIterator i, end;
  tie(i, end) = vertices(g);
  VertexSet cands(i, end);
  VertexSet nots;
  Clique clique;
  extend_clique(g, clique, cands, nots, visitor, 1);
}
