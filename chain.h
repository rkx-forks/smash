#ifndef CHAIN_H
#define CHAIN_H
#include <vector>
#include "simplex.h"

bool simplex_order(const Simplex &s1, const Simplex &s2);

class Chain {
protected:
  // NOTE : We may be assuming that all of the simplices have the same dimension.
  vector<Simplex> simplices;

public:
  Chain() : simplices() {}
  Chain(const Simplex& sx) : simplices() {
    add_simplex(sx);
  }

  template <typename Iterator>
  Chain(Iterator begin, Iterator end) : simplices(begin, end) {
    std::sort(simplices.begin(), simplices.end(), simplex_order);
  }

  bool operator==(Chain const& other);
  bool operator!=(Chain const& other);

  Chain operator+(Chain const& other) const;
  Chain operator-(Chain const& other) const;

  bool operator<(const Chain& other) const;


  vector<Simplex>::const_iterator begin() const {
    return simplices.begin();
  }
  vector<Simplex>::const_iterator end() const {
    return simplices.end();
  }

  void add_simplex(Simplex const & sx);
  bool has_simplex(Simplex const & sx);

  vector<Simplex> get_simplices();

  bool is_nonzero();
  Simplex leading_simplex();

  friend ostream& operator<<(ostream& o, Chain& c) {
    for (auto sx = c.begin(); sx != c.end(); ++sx) {
      o << *sx;
      if (sx + 1 != c.end()) 
        o << ", ";
    }
    return o;
  }
};

#endif
