#ifndef METRIC_SPACE_H
#define METRIC_SPACE_H

#include <iostream>
#include <vector>
#include <map>

#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "clique.h"

typedef float coord_t;
typedef float distance_t;

//typedef boost::multi_array<distance_t, 2> distance_matrix;
//typedef array_type::index dm_index;

typedef std::vector<coord_t> point;
typedef std::vector< std::vector<coord_t> > point_list;

using namespace boost;

class FiniteMetricSpace {
 protected:
  point_list points;
  size_t n_points;
  distance_t *distances;

  size_t dist_index(size_t i, size_t j) const {
    return i*n_points + j;
  }

  distance_t euc_dist_sq(const point &v1, const point &v2) const {
    distance_t a = 0.0;
    for (point::size_type i = 0; i < v1.size(); ++i) {
      a += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    return a;
  }

  void allocate_and_zero(point::size_type n_points) {
    distances = new distance_t[n_points*n_points];
    fill(distances, distances + n_points*n_points, 0);
  }

  FiniteMetricSpace& operator=(FiniteMetricSpace fms);
  FiniteMetricSpace(FiniteMetricSpace& fms) {}

 public:
  FiniteMetricSpace(point_list &points) :
    points(points), n_points(points.size()) {
    //    allocate_and_zero(n_points);

    /*for (point::size_type i = 0; i < n_points; ++i) 
      for (point::size_type j = 0; j < i; ++j)
        distances[dist_index(i, j)] =
        distances[dist_index(j, i)] = euc_dist_sq(points[i], points[j]);*/
  }

  distance_t get_distance(size_t i, size_t j) const {
    return euc_dist_sq(points[i], points[j]);
    //return distances[dist_index(i, j)];
  }

  size_t size() const {
    return n_points;
  }

  /* XXX todo: use filtered_graph, only create nbd graph once
   *
  NeighborhoodGraphPtr get_neighborhood_graph(coord_t epsilon_sq) {
    NeighborhoodGraphPtr ngp(new NeighborhoodGraph(size()));
    for (size_t i = 0; i < size(); ++i) 
      ngp->add_vertex(i);
    for (size_t i = 0; i < size(); ++i) {
      for (size_t j = i + 1; j < size(); ++j) {
        if (get_distance(i, j) < epsilon_sq)
          ngp->add_edge(i, j);
      }
    }
    return ngp;
    } */
};

/*****************************************/

struct NBGVertexProperties {
  int vertex;
};

struct NBGEdgeProperties {
  distance_t weight; // distance squared
};

class WeightedNeighborhoodGraph {
 protected:
  typedef adjacency_list<vecS, vecS, undirectedS,
                         NBGVertexProperties,
                         NBGEdgeProperties> nbd_graph_t;
  typedef graph_traits<nbd_graph_t>::vertex_descriptor vertex_t;
  typedef graph_traits<nbd_graph_t>::edge_descriptor edge_t;

  nbd_graph_t g;
  std::map<int, vertex_t> vertex_map;
  vector<Simplex::vertex_t> vertices;

  struct edge_pred_t { // for filtered_graph below
    edge_pred_t() : g(0) {}
    edge_pred_t(const nbd_graph_t& g, distance_t max_weight)
      : g(&g), max_weight(max_weight) {}
    bool operator()(const edge_t& edge_id) const {
      return (*g)[edge_id].weight < max_weight;
    }

    const nbd_graph_t* g;
    distance_t max_weight;
  };

  void initialize_sequential_vertices(int size) {
    for (int i = 0; i < size; ++i) {
      auto v_id = vertex(i, g);
      vertex_map[i] = v_id;
      g[v_id].vertex = i;
      vertices.push_back(i);
    }
  }
 public:
  WeightedNeighborhoodGraph(int size)
    : g(size) {
    initialize_sequential_vertices(size);
  }

  WeightedNeighborhoodGraph(const FiniteMetricSpace& fms,
                            distance_t max_epsilon)
    : g(fms.size()) {
    initialize_sequential_vertices(fms.size());
    auto max_epsilon_squared = max_epsilon * max_epsilon;

    for (size_t i = 0; i < fms.size(); ++i) {
      for (size_t j = i + 1; j < fms.size(); ++j) {
        auto d = fms.get_distance(i, j);
        bool added;
        edge_t e;

        if (d < max_epsilon_squared) {
          tie(e, added) = boost::add_edge(vertex_map[i], vertex_map[j], g);
          g[e].weight = d;
        }
      }
    }
  }

  /* construct a weighted neighborhood graph from
     a collection of simplices; all weights set to 0.0 */
  template<typename Iterator>
  WeightedNeighborhoodGraph(Iterator begin, Iterator end) {
    for (auto sx = begin; sx != end; ++sx) {
      for (auto v = sx->begin(); v != sx->end(); ++v) { // add all new vertices
        if (vertex_map.find(*v) == vertex_map.end()) { // add vertex
          auto v_id = add_vertex(g);
          vertex_map[*v] = v_id;
          g[v_id].vertex = *v;
          vertices.push_back(*v);
        }
      }
      sort(vertices.begin(), vertices.end());
      for (auto v1 = sx->begin(); v1 != sx->end(); ++v1) {
        for (auto v2 = v1 + 1; v2 != sx->end(); ++v2) {
          add_edge(*v1, *v2, 0.0);
        }
      }
    }
  }

  bool has_vertex(int i) const {
    return vertex_map.find(i) != vertex_map.end();
  }

  bool has_edge(int v1, int v2) const {
    bool found;
    edge_t e;
    auto v1_id = vertex_map.find(v1);
    auto v2_id = vertex_map.find(v2);
    assert(v1_id != vertex_map.end() &&
           v2_id != vertex_map.end());
    tie(e, found) = edge(v1_id->second, v2_id->second, g);
    return found;
  }

  distance_t get_edge_weight(int v1, int v2) const {
    bool found;
    edge_t e;
    auto v1_id = vertex_map.find(v1);
    auto v2_id = vertex_map.find(v2);
    assert(v1_id != vertex_map.end() &&
           v2_id != vertex_map.end());

    tie(e, found) = edge(v1_id->second, v2_id->second, g);
    return g[e].weight;
  }

  void add_edge(int v1, int v2, distance_t weight) {
    bool added;
    edge_t e;
    auto v1_id = vertex_map.find(v1);
    auto v2_id = vertex_map.find(v2);
    assert(v1_id != vertex_map.end() &&
           v2_id != vertex_map.end());

    tie(e, added) = boost::add_edge(v1_id->second, v2_id->second, g);
    g[e].weight = weight;
  }

  Simplex::vertex_iterator_t vertices_begin() const {
    return vertices.begin();
  }

  Simplex::vertex_iterator_t vertices_end() const {
    return vertices.end();
  }

  /* TODO: get rid of */
  template <typename CliqueFactory>
  CliqueGraphPtr get_clique_graph(CliqueFactory clique_factory,
                                  distance_t epsilon) const {
    auto eps_sq = epsilon * epsilon;
    filtered_graph<nbd_graph_t, edge_pred_t> fg(g, edge_pred_t(g, eps_sq));
    return clique_factory.make_clique_graph(fg);
  }

  friend ostream& operator<<(ostream& os, WeightedNeighborhoodGraph& ng) {
    WeightedNeighborhoodGraph::nbd_graph_t::edge_iterator e, e_end;
    os << "WeightedNeighborhoodGraph" << endl;
    os << num_vertices(ng.g) << " vertices;" << endl;
    tie(e, e_end) = edges(ng.g);
    int n_edges = 0;
    for (; e != e_end; ++e) {
      n_edges++;
      std::cout << ng.g[source(*e, ng.g)].vertex << " - "
                << ng.g[target(*e, ng.g)].vertex 
                << " (" << ng.g[*e].weight << ")" << endl;
    }
    std::cout << n_edges << " edges" << endl;
    return os;
  }
};

/*
// graph whose vertices are integers
class NeighborhoodGraph {
protected:
  typedef adjacency_list<vecS, vecS, undirectedS, NBGVertexProperties> nbd_graph_t;
  typedef graph_traits<nbd_graph_t>::vertex_descriptor vertex_t;

  nbd_graph_t g;
  std::map<int, vertex_t> vertex_map;

public:
  NeighborhoodGraph(int size)  {}
  void add_vertex(int index) { 
    vertex_t u = boost::add_vertex(g);
    vertex_map[index] = u;
    g[u].vertex = index;
  }

  void add_edge(int i1, int i2) {
    boost::add_edge(vertex_map[i1], vertex_map[i2], g);
  }

  bool has_vertex(int i) {
    return vertex_map.count(i);
  }

  bool has_edge(int i1, int i2) {
    bool found;
    graph_traits<nbd_graph_t>::edge_descriptor e;
    tie(e, found) = edge(vertex_map[i1], vertex_map[i2], g);
    return found;
  }

  template <typename CliqueFactory>
  CliqueGraphPtr get_clique_graph(CliqueFactory clique_factory) {
    return clique_factory.make_clique_graph(g);
  }
  
  void print() {
    nbd_graph_t::edge_iterator e, e_end;
    std::cout << vertex_map.size() << " vertices;" << std::endl;
    tie(e, e_end) = edges(g);
    for (; e != e_end; ++e) {
      std::cout << g[source(*e, g)].vertex << " - "
                << g[target(*e, g)].vertex << std::endl;
    }
  }
};
*/

//typedef std::shared_ptr<NeighborhoodGraph> NeighborhoodGraphPtr;


#endif
