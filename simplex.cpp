#include <assert.h>
#include <boost/functional/hash.hpp>

#include "simplex.h"


/*
void Simplex::add_vertex(vertex_t vertex) {
  assert(_dimension < MAX_VERTICES - 1);
  for (int i = _dimension + 1; i >= 0; --i) {
    if (i > 0 && vertex < vertices[i - 1]) {
      vertices[i] = vertices[i-1];
    } else {
      vertices[i] = vertex;
      break;
    }
  }
  _dimension++;
}

bool Simplex::operator==(const Simplex &other) const {
  if (_dimension != other._dimension)
    return false;
  for (int i = 0; i < _dimension + 1; ++i) {
    if (vertices[i] != other.vertices[i])
      return false;
  }
  return true;
}*/

std::ostream &operator<<(std::ostream &os, const Simplex& s) {
  os << "Simplex: ";
  for (auto i = s.begin(); i != s.end(); ++i) {
    os << *i << " ";
  }
  return os;
}

void print_simplices(const vector<Simplex>& sxs) {
  for (auto i = sxs.begin(); i != sxs.end(); ++i) {
    cout << *i << endl;
  }
}

std::size_t hash_value(Simplex const& sx) {
    size_t seed = 0;
    for (auto i = sx.begin(); i != sx.end(); ++i)
      boost::hash_combine(seed, *i);
    return seed;
}

Simplex make_simplex(vector<int> verts) {
  return Simplex(verts.begin(), verts.end());
}

/*
SubsimplexIterator subsimplices_begin(const Simplex& sx, size_t dim) {
  return SubsimplexIterator(sx, dim, false);
}

SubsimplexIterator subsimplices_end(const Simplex& sx, size_t dim) {
  return SubsimplexIterator(sx, dim, true);
  }*/
