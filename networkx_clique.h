#ifndef NETWORKX_CLIQUE_H
#define NETWORKX_CLIQUE_H

#include "clique_factory.h"
#include <boost/numeric/conversion/cast.hpp>
#include <stack>

using namespace std;

template <typename T>
bool is_sorted(T& arr) {
  return adjacent_find(arr.begin(), arr.end()) == arr.end();
}

struct NetworkxCliqueFactory {
  template <typename Graph>
  CliqueGraphPtr make_clique_graph(Graph& graph) {
    CliqueVisitor visitor;
    networkx_cliques(graph, visitor);
    return make_shared<CliqueGraph>(visitor.simplices.begin(),
                                    visitor.simplices.end());
  }
};

template <typename T>
void print_collection(string s, T& c) {
  cout << s << endl;
  for (auto i = c.begin(); i != c.end(); ++i)
    cout << *i << endl;
}

template <typename Graph, typename Visitor>
void networkx_cliques(Graph& g, Visitor& visitor) {
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertex_iterator VertexIterator;
  typedef typename graph_traits<Graph>::adjacency_iterator AdjacencyIterator;
  typedef typename graph_traits<Graph>::degree_size_type Degree;
  typedef std::vector<Vertex> VertexSet;
  typedef std::vector<Vertex> Clique;

  Degree maxconn = -1;
  map<Vertex, VertexSet> nnbrs;
  VertexSet pivotnbrs;
  Vertex max_conn_vertex;

  // create neighbor map
  VertexIterator i, i_end;
  for (tie(i, i_end) = vertices(g); i != i_end; ++i) {
    AdjacencyIterator j, j_end;
    tie(j, j_end) = adjacent_vertices(*i, g);
    VertexSet nbrs(j, j_end);
    sort(nbrs.begin(), nbrs.end());
    nnbrs[*i] = nbrs;
    if (nbrs.size() > maxconn) {
      maxconn = nbrs.size();
      max_conn_vertex = *i;
    }
  }
  pivotnbrs = nnbrs[max_conn_vertex];

  // initial set up
  tie(i, i_end) = vertices(g);
  VertexSet cand(i, i_end);
  sort(cand.begin(), cand.end());
  VertexSet smallcand;
  smallcand.reserve(cand.size());
  set_difference(cand.begin(), cand.end(),
                 pivotnbrs.begin(), pivotnbrs.end(),
                 inserter(smallcand, smallcand.begin()));
  VertexSet done;

  typedef stack< list<VertexSet> > VertexSetStack;
  vector<VertexSet> cand_stack;
  vector<VertexSet> done_stack;
  vector<VertexSet> smallcand_stack;
  Vertex n;

  Clique clique_so_far;

  while (!smallcand.empty() or !cand_stack.empty()) {
    if (!smallcand.empty()) {
      n = smallcand.back();
      smallcand.pop_back();
    } else {
      cand = cand_stack.back();
      cand_stack.pop_back();
      done = done_stack.back();
      done_stack.pop_back();
      smallcand = smallcand_stack.back();
      smallcand_stack.pop_back();
      clique_so_far.pop_back();
      continue;
    }

    clique_so_far.push_back(n);
    cand.erase(remove(cand.begin(), cand.end(), n), cand.end());
    // insert n keeping 'done' sorted
    auto ins = lower_bound(done.begin(), done.end(), n);
    done.insert(ins, n);

    VertexSet new_cand, new_done;
    auto nn_size = nnbrs[n].size();
    new_cand.reserve(min(cand.size(), nn_size));
    new_done.reserve(min(done.size(), nn_size));
    set_intersection(cand.begin(), cand.end(), 
                     nnbrs[n].begin(), nnbrs[n].end(),
                     inserter(new_cand, new_cand.begin()));
    set_intersection(done.begin(), done.end(),
                     nnbrs[n].begin(), nnbrs[n].end(),
                     inserter(new_done, new_done.begin()));
    if (new_cand.empty()) {
      if (new_done.empty())
        visitor.clique(clique_so_far, g);
      clique_so_far.pop_back();
      continue;
    }
    if (new_done.empty() && new_cand.size() == 1) {
      clique_so_far.push_back(new_cand[0]);
      visitor.clique(clique_so_far, g);
      clique_so_far.pop_back();
      clique_so_far.pop_back();
      continue;
    }
    auto numb_cand = numeric_cast<int>(new_cand.size());
    int maxconndone = -1;
    VertexSet pivotdonenbrs;
    for (auto nd = new_done.begin(); nd != new_done.end(); ++nd) {
      VertexSet cn;
      cn.reserve(min(new_cand.size(), nnbrs[*nd].size()));
      set_intersection(new_cand.begin(), new_cand.end(),
                    nnbrs[*nd].begin(), nnbrs[*nd].end(),
                    inserter(cn, cn.begin()));
      auto conn = numeric_cast<int>(cn.size());
      if (conn > maxconndone) {
        pivotdonenbrs = cn;
        maxconndone = conn;
        if (maxconndone == numb_cand)
          break;
      }
    }

    if (maxconndone == numb_cand) {
      clique_so_far.pop_back();
      continue;
    }

    int maxconn = -1;
    for (auto nc = new_cand.begin(); nc != new_cand.end(); ++nc) {
      VertexSet cn;
      cn.reserve(min(new_cand.size(), nnbrs[*nc].size()));
      set_intersection(new_cand.begin(), new_cand.end(),
                    nnbrs[*nc].begin(), nnbrs[*nc].end(),
                    inserter(cn, cn.begin()));
      auto conn = numeric_cast<int>(cn.size());
      if (conn > maxconn) {
        pivotnbrs = cn;
        maxconn = conn;
        if (maxconn == numb_cand - 1)
          break;
      }
    }

    if (maxconndone > maxconn)
      pivotnbrs = pivotdonenbrs;

    cand_stack.push_back(cand);
    done_stack.push_back(done);
    smallcand_stack.push_back(smallcand);

    cand = new_cand;
    done = new_done;

    assert(is_sorted(cand));
    assert(is_sorted(pivotnbrs));
    smallcand.clear();
    set_difference(cand.begin(), cand.end(),
                   pivotnbrs.begin(), pivotnbrs.end(),
                   inserter(smallcand, smallcand.begin()));
  }
}

template <typename Graph, typename Visitor>
void networkx_cliques_recursive(Graph& g, Visitor& visitor) {
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::vertex_iterator VertexIterator;
  typedef typename graph_traits<Graph>::adjacency_iterator AdjacencyIterator;
  typedef std::vector<Vertex> VertexSet;
  typedef std::vector<Vertex> Clique;

  map<Vertex, VertexSet> nnbrs;

  // create neighbor map
  VertexIterator i, i_end;
  for (tie(i, i_end) = vertices(g); i != i_end; ++i) {
    AdjacencyIterator j, j_end;
    tie(j, j_end) = adjacent_vertices(*i, g);
    VertexSet nbrs(j, j_end);
    sort(nbrs.begin(), nbrs.end());
    nnbrs[*i] = nbrs;
  }

  tie(i, i_end) = vertices(g);
  VertexSet cand(i, i_end);
  sort(cand.begin(), cand.end());
  VertexSet done;
  Clique clique_so_far;
  networkx_clique_extend(g, visitor, cand, done, 
                         nnbrs, clique_so_far);
}

template <typename Graph, typename Visitor, typename Container,
          typename Map, typename Clique>
void networkx_clique_extend(Graph& g, Visitor& visitor,
                            Container& cand, Container& done,
                            Map& nnbrs, 
                            Clique& clique_so_far) {
  //typedef typename graph_traits<Graph>::degree_size_type Degree;
  Container pivotnbrs;
  int maxconn = -1;
  auto n_cand = numeric_cast<int>(cand.size());
  for (auto n = done.begin(); n != done.end(); ++n) {
    Container cn;
    cn.reserve(min(cand.size(), nnbrs[*n].size()));
    set_intersection(cand.begin(), cand.end(),
                     nnbrs[*n].begin(), nnbrs[*n].end(),
                     inserter(cn, cn.begin()));
    auto conn = numeric_cast<int>(cn.size());
    if (conn > maxconn) {
      pivotnbrs = cn;
      maxconn = conn;
      if (conn == n_cand)
        return;
    }
  }
  for (auto n = cand.begin(); n != cand.end(); ++n) {
    Container cn;
    cn.reserve(min(cand.size(), nnbrs[*n].size()));
    set_intersection(cand.begin(), cand.end(),
                     nnbrs[*n].begin(), nnbrs[*n].end(),
                     inserter(cn, cn.begin()));
    auto conn = numeric_cast<int>(cn.size());
    if (conn > maxconn) {
      pivotnbrs = cn;
      maxconn = conn;
    }
  }
  Container smallercand;
  smallercand.reserve(cand.size());
  set_difference(cand.begin(), cand.end(),
                 pivotnbrs.begin(), pivotnbrs.end(),
                 inserter(smallercand, smallercand.begin()));
  for (auto n = smallercand.begin(); n != smallercand.end(); ++n) {
    cand.erase(remove(cand.begin(), cand.end(), *n), cand.end());
    /*    auto ins = lower_bound(clique_so_far.begin(),
                           clique_so_far.end(), n);
                           clique_so_far.insert(ins, n);*/
    clique_so_far.push_back(*n);
    Container new_cand, new_done;
    new_cand.reserve(min(cand.size(), nnbrs[*n].size()));
    new_done.reserve(min(done.size(), nnbrs[*n].size()));
    set_intersection(cand.begin(), cand.end(),
                     nnbrs[*n].begin(), nnbrs[*n].end(),
                     inserter(new_cand, new_cand.begin()));
    set_intersection(done.begin(), done.end(),
                     nnbrs[*n].begin(), nnbrs[*n].end(),
                     inserter(new_done, new_done.begin()));
    if (new_cand.empty() and new_done.empty()) {
      visitor.clique(clique_so_far, g);
    } else if (new_done.empty() && new_cand.size() == 1) {
      clique_so_far.push_back(new_cand[0]);
      visitor.clique(clique_so_far, g);
      clique_so_far.pop_back();
    } else {
      networkx_clique_extend(g, visitor, new_cand, 
                             new_done, nnbrs, clique_so_far);
    }
    auto last = clique_so_far.back();
    auto ins = lower_bound(done.begin(), done.end(), last);
    done.insert(ins, last);
    clique_so_far.pop_back();
  }
}

#endif
