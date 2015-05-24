#include "chain.h"

#include <algorithm>
#include <iostream>
#include <assert.h>

bool simplex_order(const Simplex &s1, const Simplex &s2) {
  return s2 < s1;
}

bool Chain::operator==(Chain const & other) {
  return (simplices == other.simplices);
}

bool Chain::operator!=(Chain const & other) {
    return !(*this == other);
}

Chain Chain::operator+(Chain const & other) const {
    return *this - other;
}

Chain Chain::operator-(Chain const & other) const {
    /* This method gets called a LOT!  We should avoid temporaries
     * and so forth.  Just construct the new chain, suggest a size
     * for its simplex vector, and call symmetric_difference.  Note that
     * the result will be sorted, so we're good on that.  RVO should 
     * eliminate any more bad copies...
     */
    Chain result;
    result.simplices.reserve(simplices.size() + other.simplices.size());

    set_symmetric_difference(begin(), end(),
                             other.begin(), other.end(),
                             back_inserter(result.simplices),
                             simplex_order);
    return result;
}

bool Chain::operator<(const Chain& other) const {
  return *simplices.begin() < *other.simplices.begin();
}

// Insert a Simplex to the chain.
// Maintain descending order of simplices.
void Chain::add_simplex(const Simplex & sx) {
  assert(!has_simplex(sx));

  auto ins = lower_bound(simplices.begin(), simplices.end(), sx, 
                         simplex_order);
  simplices.insert(ins, sx);
}

bool Chain::has_simplex(Simplex const & sx) {
  return find(simplices.begin(), simplices.end(), sx) != simplices.end();
}

vector<Simplex> Chain::get_simplices() {
    return simplices;
}

bool Chain::is_nonzero() {
    return simplices.begin() != simplices.end();
}

Simplex Chain::leading_simplex() {
    assert(is_nonzero());
    return *simplices.begin();
}
