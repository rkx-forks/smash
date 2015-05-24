#ifndef SIMPLICIAL_SET_H
#define SIMPLICIAL_SET_H

#include <iostream>
#include <unordered_map>
#include "boost/shared_ptr.hpp"
#include "boost/cast.hpp"

#include "incremental_vr2.h"
#include "chain.h"
#include "simplex.h"
#include "mytimer.h"

//using namespace boost;
using namespace std;

template <typename K, typename V>
void print_map(map<K, V> m) {
  for (auto kv = m.begin(); kv != m.end(); ++kv) {
    cout << "key=" << kv->first << ", value=" << kv->second << endl;
  }
}

class DegenerateSimplices {
 protected:
  SimplexList simplices;

  std::unordered_map<Simplex::vertex_t, Simplex::vertex_t> vertex_map;
  vector<Simplex> vertex_identifications;

  /* slower!
    void compute_identifications() {
    CliqueGraph cg(simplices.begin(), simplices.end());
    cg.connected_component_vertices(back_inserter(vertex_identifications));
    }*/

  void compute_vertex_map() {
    for (auto sx = vertex_identifications.begin();
         sx != vertex_identifications.end(); ++sx) {
      //      cout << *sx << endl;
      auto min_rep = *sx->begin();
      for (auto v = sx->begin(); v != sx->end(); ++v) {
        vertex_map[*v] = min_rep;
      }
    }
  }
  /*
     Store vertex identifications as a list of simplices
     which are pairwise non intersecting
   */
  void recompute_identifications(const Simplex &new_sx) {
    vector<Simplex> new_identifications;
    set<Simplex::vertex_t> new_verts(new_sx.begin(), new_sx.end());

    for (auto vi = vertex_identifications.begin();
         vi != vertex_identifications.end(); ++vi) {
      if (new_sx.has_intersection(*vi)) {
        copy(vi->begin(), vi->end(), inserter(new_verts, new_verts.begin()));
      } else {
        new_identifications.push_back(*vi);
      }
    }

    new_identifications.push_back(Simplex(new_verts.begin(), new_verts.end()));
    vertex_identifications = new_identifications;
  }

  DegenerateSimplices(const DegenerateSimplices&) {}
public:
  DegenerateSimplices() {}

  template <typename Iterator>
    DegenerateSimplices(Iterator start, Iterator end)
      : simplices(start, end), vertex_identifications() {
    for_each(start, end,
             [&](const Simplex &sx) { recompute_identifications(sx); });
    compute_vertex_map();
  }

  void add_simplex(Simplex& s) {
    //auto ins = lower_bound(simplices.begin(), simplices.end(), s);
    //    simplices.insert(ins, s);
    simplices.add(s);
    recompute_identifications(s);
    compute_vertex_map();
  }

  bool is_collapsed(const Simplex& s) const {
    if (s.dimension() == 0) {
      return false;
    } else {
      return simplices.contains_subsimplex(s);
    }
  }

  Simplex::vertex_t vertex_representative(Simplex::vertex_t v) const {
    auto it = vertex_map.find(v);
    if (it != vertex_map.end())
      return it->second;
    else
      return v;
  }

  template <typename Iterator>
  void all_vertex_representatives(Iterator out) const {
    set<Simplex::vertex_t> vertices;
    for (auto v = vertex_map.begin(); v != vertex_map.end(); ++v)
      vertices.insert(v->second);
    copy(vertices.begin(), vertices.end(), out);
  }

  template<typename Iterator>
  void nondegenerate_subsimplices_of(int dimension,
                                     const Simplex& s,
                                     Iterator out) const {
    if (dimension == 0) {
      set<Simplex::vertex_t> seen;
      for (auto v = s.begin(); v != s.end(); ++v) {
        auto vertex = vertex_representative(*v);
        if (seen.find(vertex) == seen.end()) {
          *out++ = Simplex(vertex);
          seen.insert(vertex);
        }
      }
    } else {
      s.each_subsimplex(dimension,
                        [&](const Simplex& ss) {
                          if (!is_collapsed(ss))
                            *out++ = ss;
                        });
    }
  }

  Chain boundary_chain(const Simplex &s) const {
    Chain c;
    if (s.dimension() == 0) {
      return c;
    } else if (s.dimension() == 1) {
      auto v1 = *(s.begin()),
        v2 = *(s.begin() + 1);
      if (vertex_representative(v1) != vertex_representative(v2)) {
        c.add_simplex(Simplex(vertex_representative(v1)));
        c.add_simplex(Simplex(vertex_representative(v2)));
      }
    } else { // dimension > 1
      s.each_subsimplex(s.dimension() - 1,
                        [&](const Simplex& ssx) {
                          if (!is_collapsed(ssx))
                            c.add_simplex(ssx);
                        });
    }
    return c;
  }

  Chain boundary_chain(const Chain &c) const {
    Chain result;
    for (auto sx = c.begin(); sx != c.end(); ++sx)
      result = result + boundary_chain(*sx);
    return result;
  }

  Chain reduce(const Chain& c) const {
    Chain result;
    for (auto sx = c.begin(); sx != c.end(); ++sx) {
      if (sx->dimension() > 0 && !is_collapsed(*sx))
        result = result + *sx;
      else if (sx->dimension() == 0) {
        auto v = *sx->begin();
        result = result + Simplex(vertex_representative(v));
      }
    }
    return result;
  }

  /*
    return true if simplex s (which is a maximal simplex)
    is acyclic in the context of this object
  */
  bool is_acyclic(const Simplex &s) const {
    if (is_collapsed(s)) {
      return true;
    } else if (s.dimension() == 0) {
      return false; // should be true?
    } else if (s.dimension() == 1) {
      auto v = s.begin();
      return vertex_representative(*v) != vertex_representative(*(++v));
    }

    int betti_current = 0,
      betti_previous = 0;
    std::unordered_map<Simplex, Chain, SimplexHasher> pivots;
    for (int d = 0; d < s.dimension() + 1; ++d) {
      vector<Simplex> subsimplices; // nondegenerate subsimplices of s
      nondegenerate_subsimplices_of(d, s, back_inserter(subsimplices));
      for (auto subsx = subsimplices.begin();
           subsx != subsimplices.end(); ++subsx) {
        auto boundary = boundary_chain(*subsx);
        while (boundary.is_nonzero() &&
               pivots.count(boundary.leading_simplex()))
          boundary = boundary - pivots[boundary.leading_simplex()];
        if (boundary.is_nonzero()) {
          pivots[boundary.leading_simplex()] = boundary;
          if (d > 0)
            betti_previous -= 1;
        } else {
          betti_current += 1;
        }
      }

      if (d > 1 && betti_previous != 0) {
        return false;
      }
      betti_previous = betti_current;
      betti_current = 0;
    }

    return betti_previous == 0;
  }

  void print() const { simplices.print(); }

  void print_info() const {
    cout << "degenerate simplices: " << endl;
    simplices.print_info();
  }
};

/************************************/

class SimplicialSet {
protected:
  vector<Simplex> maximal_simplices;
  const DegenerateSimplices collapsed_simplices;

  // populated in get_homology
  vector< vector<Chain> > homology_basis;

  typedef vector<Chain>::iterator chain_iterator_t;

  template <typename Iterator>
  void nondegenerate_subsimplices_dim(int dim, Iterator out) const {
    set<Simplex> sxs;
    for (auto sx = maximal_simplices.begin();
         sx != maximal_simplices.end(); ++sx) {
      if (dim <= sx->dimension())
        collapsed_simplices.nondegenerate_subsimplices_of
          (dim, *sx, inserter(sxs, sxs.begin()));
    }

    if (dim == 0) { // add ALL vertices, even if they belong to collasped sxs
      vector<Simplex::vertex_t> verts;
      collapsed_simplices.all_vertex_representatives(back_inserter(verts));
      for_each(verts.begin(), verts.end(),
               [&](Simplex::vertex_t v) { sxs.insert(Simplex(v)); });
    }

    copy(sxs.begin(), sxs.end(), out);
  }

  template <typename Iterator>
  void old_nondegenerate_subsimplices(int max_dim, Iterator out) const {
    for (int i = 0; i < max_dim + 1; ++i)
      nondegenerate_subsimplices_dim(i, out);
  }

  template <typename Iterator>
  void nondegenerate_subsimplices(int max_dim, Iterator out) const {
    WeightedNeighborhoodGraph ng(maximal_simplices.begin(),
                                 maximal_simplices.end());
    // get 0-simplices
    set<Simplex> zero_sxs;
    for (auto s = maximal_simplices.begin();
         s != maximal_simplices.end(); ++s) {
      collapsed_simplices.nondegenerate_subsimplices_of
          (0, *s, inserter(zero_sxs, zero_sxs.begin()));
    }
    vector<Simplex::vertex_t> verts;
    collapsed_simplices.all_vertex_representatives
      (back_inserter(verts));
    for_each(verts.begin(), verts.end(),
             [&](Simplex::vertex_t v)
             { zero_sxs.insert(Simplex(v)); });

    copy(zero_sxs.begin(), zero_sxs.end(), out);

    // get all n-simplices, n>0
    auto visitor = make_simplicial_vr_visitor(collapsed_simplices);
    incremental_vr(ng, max_dim, visitor);
    visitor.all_simplices(out);
  }

public:
  SimplicialSet(vector<Simplex> maximal_simplices,
                vector<Simplex> collapsed_simplices) :
    maximal_simplices(maximal_simplices),
    collapsed_simplices(collapsed_simplices.begin(),
                        collapsed_simplices.end()),
    homology_basis() {}

  bool is_collapsed(const Simplex& s) const {
    return collapsed_simplices.is_collapsed(s);
  }

  Chain reduce(const Chain& c) const {
    return collapsed_simplices.reduce(c);
  }

  Simplex::vertex_t vertex_representative(Simplex::vertex_t v) const {
    return collapsed_simplices.vertex_representative(v);
  }

  void print_hom_pivots(vector< map<Simplex, Chain> > hpivots) const {
    for (size_t i = 0; i < hpivots.size(); ++i) {
      cout << "dim " << i << endl;
      auto pivs = hpivots[i];
      for (auto it = pivs.begin(); it != pivs.end(); ++it) {
        cout << it->first << " ----> ";
        cout << it->second << endl;
      }
    }
  }

  void compute_homology(int max_dim) {
      std::unordered_map<Simplex, Chain, SimplexHasher> image_pivots;
      std::unordered_map<Simplex, Chain, SimplexHasher> kernel_pivots;
      std::unordered_map<Simplex, Chain, SimplexHasher> cokernel_pivots;
      vector<int> betti(max_dim + 1);
      vector< std::unordered_map<Simplex, Chain, SimplexHasher> >
          homology_pivots(max_dim + 1);

    vector<Simplex> nondegenerate_sxs;

    TIMER("nondeg.")
    nondegenerate_subsimplices(max_dim,
                               back_inserter(nondegenerate_sxs));
    ENDTIMER;

    /*cout << "mamximal: ";
    print_simplices(maximal_simplices);
    cout << "collapsed: ";
    collapsed_simplices.print();
    cout << "nondegenerate: ";
    print_simplices(nondegenerate_sxs);
    cout << "XXX" << endl;*/

    for (auto sx = nondegenerate_sxs.begin();
         sx != nondegenerate_sxs.end(); ++sx) {
      auto boundary = collapsed_simplices.boundary_chain(*sx);
      Chain preimage(*sx);
      int d = sx->dimension();
      //cout << "current sx " << *sx << endl;
      while (boundary.is_nonzero() &&
             image_pivots.count(boundary.leading_simplex())) {
        //assert(collapsed_simplices.boundary_chain(preimage) ==
        //               boundary);
        auto s = boundary.leading_simplex();
        boundary = boundary - image_pivots[s];
        preimage = preimage - cokernel_pivots[s];
        //assert(collapsed_simplices.boundary_chain(preimage) ==
        //       boundary);
      }
      /*cout << "boundary ";
      boundary.print();
      cout << "preimage ";
      preimage.print();
      cout << endl;*/

      if (boundary.is_nonzero()) {
        image_pivots[boundary.leading_simplex()] = boundary;
        cokernel_pivots[boundary.leading_simplex()] = preimage;
        if (d >= 1) {
          auto it = homology_pivots[d - 1].find(boundary.leading_simplex());
          assert(it != homology_pivots[d - 1].end());
          homology_pivots[d - 1].erase(it);
          //print_hom_pivots(homology_pivots);
        }
      } else { // boundary == 0
        homology_pivots[d][preimage.leading_simplex()] = preimage;
        //print_hom_pivots(homology_pivots);
      }
    }

    /*    homology_basis.resize(max_dim);
    for (int d = 0; d < max_dim; ++d) {
      for (auto kv = homology_pivots[d].begin();
           kv != homology_pivots[d].end(); ++kv) {
        homology_basis[d].push_back(kv->second);
      }
      sort(homology_basis[d].begin(), homology_basis[d].end());
      }*/

    homology_basis.resize(max_dim);
    for (int d = 0; d < max_dim; ++d) {
      for (auto kv = homology_pivots[d].begin();
           kv != homology_pivots[d].end(); ++kv) {
        assert((kv->second).is_nonzero());
        auto ins = lower_bound(homology_basis[d].begin(),
                               homology_basis[d].end(), kv->second);
        homology_basis[d].insert(ins, kv->second);
      }
    }
  }

  chain_iterator_t homology_basis_begin(int dim) {
    return homology_basis[dim].begin();
  }

  chain_iterator_t homology_basis_end(int dim) {
    return homology_basis[dim].end();
  }

  // provisional, rough measure of size
  int rough_size(void) const {
      unordered_set<Simplex> sxs;
      for (auto i = maximal_simplices.begin();
              i != maximal_simplices.end(); ++i) {
          vector<Simplex> tmp;
          i->all_subsimplices(3, back_inserter(tmp));
          for_each(tmp.begin(), tmp.end(),
                  [&](Simplex s) { if (!is_collapsed(s)) sxs.insert(s); });
      }
      return numeric_cast<int>(sxs.size());
  }

  /*int uncollapsed_size(void) const {
    set<Simplex> sxs;
    for (auto i = maximal_simplices.begin();
         i != maximal_simplices.end(); ++i) {
      i->all_subsimplices(3, inserter(sxs, sxs.end()));
    }
    for (auto i = collapsed_simplices.begin();
         i != collapsed_simplices.end(); ++i) {
      i->all_subsimplices(3, inserter(sxs, sxs.end()));
    }
    return numeric_cast<int>(sxs.size());
    }*/

  void print_info(void) const {
    cout << "Simplicial Set, maximal simplices:" << endl;
    cout << "max_sxs size " << maximal_simplices.size();
    cout << ", total sxs approx. " << rough_size();
    /*cout << ", uncollapsed " << uncollapsed_size() << endl;*/

    print_simplices(maximal_simplices);
    cout << "collapsed simplices:" << endl;
    collapsed_simplices.print();
    collapsed_simplices.print_info();
  }

  void print_betti(int dim) {
    cout << "b" << dim << "=" <<
      distance(homology_basis_begin(dim),
               homology_basis_end(dim)) << endl;
  }
};

typedef std::shared_ptr<SimplicialSet> SimplicialSetPtr;

class FilteredSimplicialSet {
  vector<SimplicialSetPtr> simplicial_sets;
  vector<coord_t> epsilons;
  int max_dim;

  typedef vector<SimplicialSetPtr>::size_type vector_size_t;

  // induced_maps[t][dim][c] gives the value of the induced map
  // at time t in dimension dim on the chain c
  vector< vector< map<Chain, Chain> > > induced_maps;
  // p_intervals[dim] is a list of p_intervals in dimension
  // dim using the right-open interval convention
  vector< vector< pair<int, int> > > p_intervals;

  void compute_induced_maps() {
    induced_maps.resize(filtration_length() - 1);

    for (int i = 0; i < filtration_length() - 1; ++i) {
      auto domain = simplicial_sets[i];
      auto codomain = simplicial_sets[i + 1];
      induced_maps[i].resize(max_dim);

      //cout << "induced map at level " << i << endl;

      for (int dim = 0; dim < max_dim; ++dim) {
        //cout << " dim " << dim << ":" << endl;
        for (auto e = domain->homology_basis_begin(dim);
             e != domain->homology_basis_end(dim); ++e) {
          induced_maps[i][dim][*e] = codomain->reduce(*e);
          //cout << "   " << *e << " -> " << induced_maps[i][dim][*e] << endl;
        }
      }
    }
  }

  void compute_p_intervals(int dim) {
    assert(dim <= max_dim);
    assert(simplicial_sets.size() > 0);

    std::unordered_map<Simplex, int, SimplexHasher> birth_times;
    vector<Chain> image;
    std::unordered_map<Simplex, Chain, SimplexHasher> pivots;

    // compute basis for first space
    for (auto c = simplicial_sets[0]->homology_basis_begin(dim);
         c != simplicial_sets[0]->homology_basis_end(dim); ++c) {
      auto chain = *c;
      while (chain.is_nonzero() &&
             pivots.count(chain.leading_simplex()))
        chain = chain - pivots[chain.leading_simplex()];
      if (chain.is_nonzero()){
        pivots[chain.leading_simplex()] = chain;
        birth_times[chain.leading_simplex()] = 0;
        image.push_back(chain);
      }
    }

    // main loop
    for (int i = 1; i < filtration_length(); ++i) {
      auto s = simplicial_sets[i];
      auto& f = induced_maps[i - 1];

      pivots.clear();
      std::unordered_map<Simplex, int, SimplexHasher> new_times;
      vector<Chain> new_image;

      for (auto ch = image.begin(); ch != image.end(); ++ch) {
        auto im = f[dim][*ch];
        while (im.is_nonzero() &&
               pivots.count(im.leading_simplex()))
          im = im - pivots[im.leading_simplex()];
        if (im.is_nonzero()) {
          pivots[im.leading_simplex()] = im;
          assert(birth_times.count(ch->leading_simplex()));
          new_times[im.leading_simplex()] =
            birth_times[ch->leading_simplex()];
          new_image.push_back(im);
        } else {
          auto p = make_pair(birth_times[ch->leading_simplex()], i);
          p_intervals[dim].push_back(p);
        }
      }

      for (auto ch = s->homology_basis_begin(dim);
           ch != s->homology_basis_end(dim); ++ch) {
        auto chain = *ch;
        while (chain.is_nonzero() &&
               pivots.count(chain.leading_simplex()))
          chain = chain - pivots[chain.leading_simplex()];
        if (chain.is_nonzero()) {
          pivots[chain.leading_simplex()] = chain;
          new_times[chain.leading_simplex()] = i;
          new_image.push_back(chain);
        }
      }

      birth_times = new_times;
      image = new_image;
    }

    for (auto kv = birth_times.begin();
         kv != birth_times.end(); ++kv) {
      auto p = make_pair(kv->second, -1);
      p_intervals[dim].push_back(p);
    }
  }

 public:
  FilteredSimplicialSet(int max_dim) : max_dim(max_dim) {}

  int filtration_length() const {
    return static_cast<int>(simplicial_sets.size());
  }

  void add_simplicial_set(coord_t eps, SimplicialSetPtr ss) {
    assert(!epsilons.size() || eps < epsilons.back());
    simplicial_sets.insert(simplicial_sets.begin(), ss);
    epsilons.insert(epsilons.begin(), eps);
  }

  SimplicialSetPtr last_simplical_set(void) {
    return simplicial_sets.back();
  }


  void compute_p_intervals() {
    for (auto ss = simplicial_sets.begin();
         ss != simplicial_sets.end(); ++ss) {
      (*ss)->compute_homology(max_dim);
    }
    compute_induced_maps();

    p_intervals.resize(max_dim);
    for (int dim = 0; dim < max_dim; ++dim) {
      compute_p_intervals(dim);
    }
  }


  void print_info() {
    auto e = epsilons.begin();
    auto ss = simplicial_sets.begin();
    for (;e != epsilons.end(), ss != simplicial_sets.end();
         ++e, ++ss) {
      cout << "e=" << *e << " SimplicialSet w/ roughly "
           << (*ss)->rough_size() << " sxs" << endl;
      //(*ss)->print_info();
      for (int d = 0; d < max_dim; ++d) {
        (*ss)->print_betti(d);
      }
    }
  }

  void print_induced_maps() {
    cout << "induced maps" << endl;
    for (int i = 0; i < filtration_length() - 1; ++i) {
      cout << "time " << i << endl;
      for (int d = 0; d < max_dim; ++d) {
        cout << " dim " << d << endl;
        for (auto kv = induced_maps[i][d].begin();
             kv != induced_maps[i][d].end(); ++kv) {
          Chain c1 = kv->first;
          Chain c2 = kv->second;
          cout << "  " << c1 << " -> " <<
            c2 << endl;
        }
      }
    }
  }

  void print_p_intervals() {
    for (int d = 0; d < max_dim; ++d) {
      cout << "dimension " << d << endl;
      for (auto i = p_intervals[d].begin();
           i != p_intervals[d].end(); ++i)
        cout << " (" << i->first << ", " << i->second << ")" << endl;
    }
  }
};

typedef std::shared_ptr<FilteredSimplicialSet> FilteredSimplicialSetPtr;

#endif
