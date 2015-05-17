#ifndef _INCREMENTAL_VR2_H_
#define _INCREMENTAL_VR2_H_

#include <boost/cast.hpp>

#include <vector>
#include <unordered_map>

class LowerNeighbors {
protected:
  typedef vector<int>::const_iterator neighbors_iterator_t;
  std::unordered_map<int, vector<int>> neighbors;
public:
  LowerNeighbors(WeightedNeighborhoodGraph& ng) 
    : neighbors() {
    for (auto v1 = ng.vertices_begin(); v1 != ng.vertices_end(); ++v1) {
      
      for (auto v2 = v1 + 1; v2 != ng.vertices_end(); ++v2) {
        assert(*v1 < *v2);
        if (ng.has_edge(*v1, *v2))
          neighbors[*v1].push_back(*v2);
      }
      sort(neighbors[*v1].begin(), neighbors[*v1].end());
    }
  }

  neighbors_iterator_t begin(int n) {
    return neighbors[n].begin();
  }

  neighbors_iterator_t end(int n) {
    return neighbors[n].end();
  }
};

template <typename Graph, typename Visitor>
void incremental_vr(Graph& g, int k, Visitor& visitor) {
  Simplex sx, zero;
  LowerNeighbors ln(g);

  for (auto v1 = g.vertices_begin(); v1 != g.vertices_end(); ++v1) {
    vector<int> verts;
    verts.push_back(*v1);
    add_cofaces(g, ln, k, verts, ln.begin(*v1), ln.end(*v1), 
                0.0, visitor);
  }
}

template <typename Graph, typename Collection,
          typename InIterator, typename Visitor>
  void add_cofaces(Graph& g, LowerNeighbors& ln, int k, 
                   Collection& t, InIterator N_begin, InIterator N_end,
                   distance_t weight, Visitor& visitor) {

  if (t.size() >= 1) {
    Simplex s(t.begin(), t.end());
    visitor.visit(s, weight);
  }

  if (numeric_cast<int>(t.size()) < k + 1) {
    for (auto i = N_begin; i != N_end; ++i) {
      vector<int> verts(t.begin(), t.end());
      vector<int> NcapLowerNbrs;
      distance_t weight = visitor.get_weight(g, verts, *i);
      
      verts.push_back(*i);
      assert(is_sorted(ln.begin(*i), ln.end(*i)));
      assert(is_sorted(N_begin, N_end));
      set_intersection(ln.begin(*i), ln.end(*i),
                       N_begin, N_end, 
                       back_inserter(NcapLowerNbrs));
      
      add_cofaces(g, ln, k, verts, NcapLowerNbrs.begin(),
                  NcapLowerNbrs.end(), weight, visitor);
    }
  }
}

template <typename DS>
class SimplicialVrVisitor {
    std::unordered_map<int, vector<Simplex>> sx_map;
    DS& degenerate;
    int max_dim;
public:
    SimplicialVrVisitor(DS& ds) 
        : degenerate(ds) {}

  void visit(Simplex& s, distance_t weight) {
    if (s.dimension() > 0 &&
        !degenerate.is_collapsed(s)) {
      sx_map[s.dimension()].push_back(s);
      max_dim = max(max_dim, s.dimension());
    }
  }

  template <typename Collection>
  distance_t get_weight(WeightedNeighborhoodGraph& ng, 
                  Collection& c, Simplex::vertex_t new_vert) {
    return 0.0;
  }

  template <typename Iterator>
  void all_simplices(Iterator out) {
    for (int d = 1; d <= max_dim; ++d) {
      sort(sx_map[d].begin(), sx_map[d].end());
      copy(sx_map[d].begin(), sx_map[d].end(), out);
    }
  }
};

template <typename DS>
SimplicialVrVisitor<DS> 
make_simplicial_vr_visitor(DS& degenerate) {
  return SimplicialVrVisitor<DS>(degenerate);
}

class TestingVrVisitor {
    std::unordered_map<int, vector<Simplex> > sx_map;
    std::unordered_map<Simplex, distance_t, SimplexHasher> weight_map;
  int max_dim;
  vector<Simplex> visited_order;
public:
  void visit(Simplex& s, distance_t weight) {
    max_dim = max(max_dim, s.dimension());
    sx_map[s.dimension()].push_back(s);
    weight_map[s] = weight;
    visited_order.push_back(s);
  }

  template <typename Collection>
  distance_t get_weight(WeightedNeighborhoodGraph& ng, 
                  Collection& c, Simplex::vertex_t new_vert) {
    return 0.0;
  }

  int num_sxs_of_dimension(int dim) {
    return sx_map[dim].size();
  }

  void print() {
    for (auto s = visited_order.begin();
         s != visited_order.end(); ++s) {
      cout << *s << endl;
    }
    /*    for (int d = 0; d <= max_dim; ++d) {
      sort(sx_map[d].begin(), sx_map[d].end());
      for (auto sx = sx_map[d].begin(); 
           sx != sx_map[d].end(); ++sx) {
        cout << *sx << endl;
      }
      }*/
  }
};

#endif /* _INCREMENTAL_VR2_H_ */
