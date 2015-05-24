#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <cassert>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <memory>
#include <algorithm>

#include "comb.h"


using namespace std;

class Simplex {
 public:
  typedef int vertex_t;
  typedef vector<vertex_t> vertex_vector_t;
  typedef vertex_vector_t::const_iterator vertex_iterator_t;

 protected:
  shared_ptr<vertex_vector_t> vertices;

 public:
  Simplex() : vertices(make_shared<vertex_vector_t>()) {}

  Simplex(vertex_t v) : vertices(make_shared<vertex_vector_t>()) {
    vertices->push_back(v);
  }

  template <typename Iterator>
  Simplex(Iterator begin, Iterator end) : vertices() {
    vertices = make_shared<vertex_vector_t>(begin, end);
    sort(vertices->begin(), vertices->end());
    assert(vertices->size() > 0);
  }

  Simplex(const Simplex &s) {
    vertices = s.vertices;
  }

  Simplex& operator=(Simplex const& other) {
    vertices = other.vertices;
    return *this;
  }

  vertex_iterator_t begin() const { return vertices->begin(); }
  vertex_iterator_t end() const { return vertices->end(); }

  /*  void add_vertex(vertex_t v) {
    auto ins = lower_bound(vertices.begin(), vertices.end(), v);
    vertices.insert(ins, v);
    }*/

  int dimension() const { return vertices->size() - 1; }

  bool has_vertex(vertex_t v) const {
    return binary_search(begin(), end(), v);
  }

  /*  this is subsimplex of other
   */
  bool is_subsimplex(const Simplex& other) const {
    return includes(other.begin(), other.end(),
                    begin(), end());
  }

  bool has_intersection(const Simplex& other) const {
    vertex_iterator_t i = begin(),
      j = other.begin(),
      i_end = end(),
      j_end = other.end();

    while (i != i_end && j != j_end) {
      if (*i < *j) ++i;
      else if (*j < *i) ++j;
      else {
        return true;
      }
    }
    return false;
  }

  template <typename Function>
  void each_subsimplex(size_t dim, Function f) const {
    assert(dim <= static_cast<size_t>(dimension()));
    vector<vertex_t> verts(begin(), end());
    typedef vector<vertex_t>::iterator v_it;

    auto begin = verts.begin();
    auto mid = verts.begin() + dim + 1;
    auto end = verts.end();
    for_each_combination(begin, mid, end,
                         [&](v_it b, v_it e)
                         { f(Simplex(b, e)); return false; });
  }

  template <typename Iterator>
  void all_subsimplices(size_t max_dimension, Iterator out) const {
    size_t max_dim = min(max_dimension, (size_t)dimension());
    for (size_t dim = 0; dim < max_dim + 1; ++dim) {
      each_subsimplex(dim, [&](const Simplex& s) { *out++ = s; });
    }
  }

  bool operator==(const Simplex& other) const {
    return (vertices == other.vertices ||
            *vertices == *other.vertices);
  }
  bool operator!=(Simplex& other) { return !(*this == other); }
  bool operator<(const Simplex& other) const {
    //    if (dimension() < other.dimension())
    //      return true;
    return lexicographical_compare(begin(), end(),
                                   other.begin(),
                                   other.end());
  }
  friend std::ostream &operator<<(std::ostream &os, const Simplex& s);
};

size_t hash_value(Simplex const& p);
void print_simplices(const vector<Simplex>& sxs);
Simplex make_simplex(vector<int> verts);

struct SimplexHasher {
    std::size_t operator()(const Simplex& s) const {
        using std::size_t;
        using std::hash;
        using std::string;

        std::size_t hash_val = 0;
        for (auto i = s.begin(); i != s.end(); ++i) {
            hash_val ^= hash<Simplex::vertex_t>()(*i);
        }
        return hash_val;
    }
};

//////////////////////////////////////////////
//
class SimplexList {
  unordered_map<Simplex::vertex_t, vector<Simplex> > sx_map;
  //  vector<Simplex> simplices;

  typedef vector<Simplex>::const_iterator simplex_iterator_t;

 public:
  SimplexList() {}

  template <typename Iterator>
  SimplexList(Iterator begin, Iterator end) {
    for_each(begin, end, [&](const Simplex &sx) {
        add(sx);
      });
  }

  void add(const Simplex &sx) {
    for (auto v = sx.begin(); v != sx.end(); ++v) {
      sx_map[*v].push_back(sx);
    }
    //simplices.push_back(sx);
  }

  // can get rid of?
  /*
  simplex_iterator_t begin() const { return simplices.begin(); }
  simplex_iterator_t end() const { return simplices.end(); }*/

  bool contains_subsimplex(const Simplex &sx) const {
    auto v = *sx.begin();
    if (!sx_map.count(v))
      return false;

    auto sx_list = sx_map.at(v);
    return any_of(sx_list.begin(), sx_list.end(),
                  [&](const Simplex& face)
                  { return sx.is_subsimplex(face); });
  }

  void print() const {
    cout << "(SimplexList with some simplices)" << endl;
    //print_simplices(simplices);
  }

  void print_info() const {
    cout << "simplex list:" << endl;
    for (auto kv = sx_map.begin(); kv != sx_map.end(); ++kv) {
      cout << " key " << kv->first << " len " << (kv->second).size() << endl;
    }
  }
};

#endif
