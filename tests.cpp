#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Hello

#include <boost/test/unit_test.hpp>
#include <boost/array.hpp>

#include "simplex.h"
#include "metric_space.h"
#include "incremental_vr2.h"
#include "clique_factory.h"
#include "simplicial_set.h"

BOOST_AUTO_TEST_SUITE(SimplexTests)

BOOST_AUTO_TEST_CASE(SimplexCreate) {
  Simplex s1, s2;
  BOOST_CHECK(s1 == s2);
  int s3v[] = {22};
  Simplex s3(s3v, s3v + 1);
  BOOST_CHECK(s1 != s3);
  int s4v[] = {22};
  Simplex s4(s4v, s4v + 1);
  BOOST_CHECK(s3 == s4);
}

BOOST_AUTO_TEST_CASE(SimplexIntersection) {
  Simplex zero, s1 = make_simplex({7, 11, 23}),
    s2 = make_simplex({11, 12, 23}),
    s3 = make_simplex({6, 8, 18});

  BOOST_CHECK(!zero.has_intersection(s1));
  BOOST_CHECK(s1.has_intersection(s2));
  BOOST_CHECK(s2.has_intersection(s1));
  BOOST_CHECK(!s1.has_intersection(s3));
  BOOST_CHECK(!s3.has_intersection(s1));

  BOOST_CHECK(zero < s2);
  BOOST_CHECK(zero < s2);
  BOOST_CHECK(s1 < s2);
  BOOST_CHECK(s3 < s1);
}

BOOST_AUTO_TEST_CASE(SimplexSubsimplices) {
  Simplex s = make_simplex({1, 2, 3});
  set<Simplex> sxs;

#define SXS_CHECK_SIZE(dim, expected) \
  sxs.clear();                                                      \
  s.each_subsimplex(dim, [&](const Simplex& s) { sxs.insert(s); }); \
  BOOST_CHECK_EQUAL(sxs.size(), (size_t)expected);
#define ALL_SXS_CHECK_SIZE(dim, expected)       \
  sxs.clear();                                       \
  s.all_subsimplices(dim, inserter(sxs, sxs.end())); \
  BOOST_CHECK_EQUAL(sxs.size(), (size_t)expected);
  
  
  SXS_CHECK_SIZE(0, 3);
  SXS_CHECK_SIZE(1, 3);
  SXS_CHECK_SIZE(2, 1);

  ALL_SXS_CHECK_SIZE(0, 3);
  ALL_SXS_CHECK_SIZE(1, 6);
  ALL_SXS_CHECK_SIZE(2, 7);
  ALL_SXS_CHECK_SIZE(3, 7);
}

BOOST_AUTO_TEST_CASE(SimplexIsSubsimplex) {
  Simplex s1 = make_simplex({1, 2, 3, 4, 5}),
    s2 = make_simplex({2, 4, 5});
  BOOST_CHECK(s2.is_subsimplex(s1));
}

BOOST_AUTO_TEST_CASE(SimplexComparison) {
  Simplex s1 = make_simplex({1, 2, 3, 4, 5}),
    s2 = make_simplex({2, 4, 5});
  BOOST_CHECK(s1 < s2);
}

BOOST_AUTO_TEST_SUITE_END();

struct MetricSpaceFixture {
  point_list pts;
  FiniteMetricSpace *fms;

  MetricSpaceFixture() {
    point p1(2), p2(2), p3(2), p4(2);
    p1[0] = 0.0; p1[1] = 0.0;
    p2[0] = 1.0; p2[1] = 0.0;
    p3[0] = -10.0; p3[1] = 0.0;
    p4[0] = 10.0; p4[1] = 10.0;
    
    pts.push_back(p1);
    pts.push_back(p2);
    pts.push_back(p3);  
    pts.push_back(p4);
    fms = new FiniteMetricSpace(pts);
  }
  ~MetricSpaceFixture() { delete fms; }
};

BOOST_FIXTURE_TEST_SUITE(MetricSpace, MetricSpaceFixture)

BOOST_AUTO_TEST_CASE(MetricSpaceCreate) {
  BOOST_CHECK_EQUAL(fms->size(), 4);
  BOOST_CHECK_EQUAL(fms->get_distance(0, 0), 0.0);
  BOOST_CHECK_EQUAL(fms->get_distance(0, 1), 1.0);
  BOOST_CHECK_EQUAL(fms->get_distance(1, 0), 1.0);
  BOOST_CHECK_EQUAL(fms->get_distance(0, 2), 100.0);
  BOOST_CHECK_EQUAL(fms->get_distance(2, 0), 100.0);
}

BOOST_AUTO_TEST_CASE(NeighborhoodEdgeTests) {
  WeightedNeighborhoodGraph ng(4);

  ng.add_edge(0, 1, 1.0);
  ng.add_edge(1, 2, 1.0);
  ng.add_edge(2, 3, 1.0);
  ng.add_edge(3, 0, 1.0);

  for (int i = 0; i < 7; ++i) {
    if (i < 4)
      BOOST_CHECK(ng.has_vertex(i));
    else
      BOOST_CHECK(!ng.has_vertex(i));
  }

  BOOST_CHECK(ng.has_edge(0, 1));
  BOOST_CHECK(ng.has_edge(0, 3));
  BOOST_CHECK(ng.has_edge(3, 0));
  BOOST_CHECK(!ng.has_edge(0, 2));
  BOOST_CHECK(!ng.has_edge(2, 0));
}

BOOST_AUTO_TEST_CASE(NeighborhoodSortedVertices) {
  WeightedNeighborhoodGraph ng(4);

  ng.add_edge(0, 1, 1.0);
  ng.add_edge(1, 2, 1.0);
  ng.add_edge(2, 3, 1.0);
  ng.add_edge(3, 0, 1.0);
  
  BOOST_CHECK(is_sorted(ng.vertices_begin(), ng.vertices_end()));
}

BOOST_AUTO_TEST_CASE(CliqueGraphTest) {
  WeightedNeighborhoodGraph ng(4);
  Simplex s1 = make_simplex({0, 1}),
    s2 = make_simplex({2, 3}),
    s3 = make_simplex({1, 2}),
    s4 = make_simplex({0, 3}),
    s5 = make_simplex({0, 2}),
    s6 = make_simplex({1, 3});

  ng.add_edge(0, 1, 1.0);
  ng.add_edge(1, 2, 1.0);
  ng.add_edge(3, 2, 1.0);
  ng.add_edge(3, 0, 1.0);

  CliqueGraphPtr cgp = ng.get_clique_graph(BKCliqueFactory(), 2.0);
  //cgp->print_verts();
  BOOST_CHECK(cgp->has_vertex(s1));
  BOOST_CHECK(cgp->has_vertex(s2));
  BOOST_CHECK(cgp->has_vertex(s3));
  BOOST_CHECK(cgp->has_vertex(s4));
  BOOST_CHECK(!cgp->has_vertex(s5));
  BOOST_CHECK(!cgp->has_vertex(s6));
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE(IncrementalVR)

BOOST_AUTO_TEST_CASE(LowerNeighbors) {
  WeightedNeighborhoodGraph ng(3);
  ng.add_edge(0, 1, 1.0);
  ng.add_edge(1, 2, 1.0);
  ng.add_edge(0, 2, 1.0);
  
  TestingVrVisitor visitor;
  incremental_vr(ng, 3, visitor);
  BOOST_CHECK_EQUAL(visitor.num_sxs_of_dimension(0), 3);
  BOOST_CHECK_EQUAL(visitor.num_sxs_of_dimension(1), 3);
  BOOST_CHECK_EQUAL(visitor.num_sxs_of_dimension(2), 1);
}

BOOST_AUTO_TEST_CASE(NeighborhoodFromSimplexes) {
  vector<Simplex> v;
  v.push_back(make_simplex({0, 2, 3}));
  v.push_back(make_simplex({3, 5}));

  WeightedNeighborhoodGraph ng(v.begin(), v.end());
  TestingVrVisitor vis;

  incremental_vr(ng, 3, vis);
  BOOST_CHECK_EQUAL(vis.num_sxs_of_dimension(0), 4);
  BOOST_CHECK_EQUAL(vis.num_sxs_of_dimension(1), 4);
  BOOST_CHECK_EQUAL(vis.num_sxs_of_dimension(2), 1);
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE(SimplicialSetTests)
BOOST_AUTO_TEST_CASE(VertexIdentificationTests) {
  int s1v[] = {0, 1}, s2v[] = {1, 2}, s3v[] = {2, 3};
  vector<Simplex> maximal;
  vector<Simplex> collapsed;
  collapsed.push_back(Simplex(s1v, s1v+2));
  collapsed.push_back(Simplex(s2v, s2v+2));
  collapsed.push_back(Simplex(s3v, s3v+2));

  SimplicialSet ss(maximal, collapsed);
  BOOST_CHECK_EQUAL(ss.vertex_representative(0), 0);
  BOOST_CHECK_EQUAL(ss.vertex_representative(1), 0);
  BOOST_CHECK_EQUAL(ss.vertex_representative(2), 0);
  BOOST_CHECK_EQUAL(ss.vertex_representative(3), 0);
  BOOST_CHECK_EQUAL(ss.vertex_representative(4), 4);
}

BOOST_AUTO_TEST_CASE(VertexIdentificationBug1) {
  vector<Simplex> collapsed;
  collapsed.push_back(make_simplex({1, 3}));
  collapsed.push_back(make_simplex({2, 4}));
  collapsed.push_back(make_simplex({3, 4}));

  DegenerateSimplices d(collapsed.begin(), collapsed.end());
  BOOST_CHECK_EQUAL(d.vertex_representative(2), 1);
}

BOOST_AUTO_TEST_CASE(VertexIdentificationBug2) {
  vector<Simplex> collapsed;

  collapsed.push_back(make_simplex({0, 1, 9}));
  collapsed.push_back(make_simplex({2, 5}));
  collapsed.push_back(make_simplex({6, 8, 9}));

  DegenerateSimplices d(collapsed.begin(), collapsed.end());
  BOOST_CHECK_EQUAL(d.vertex_representative(8), 0);
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE(ChainTests) 
BOOST_AUTO_TEST_CASE(ChainAdd) {
  Simplex::vertex_t points[] = {1,2,3,4,5};
  Simplex sx_list[] = {
    Simplex(points, points+2), Simplex(points+1, points+4),
    Simplex(points+3, points+5), Simplex(points+2, points+4)
  };
  
  Chain ch1(sx_list, sx_list+2);
  Chain ch2(sx_list+2, sx_list+4);
  Chain c3 = ch1 + ch2;
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE(DegenerateSimplicesTests)
BOOST_AUTO_TEST_CASE(IsCollapsedTests) {
  Simplex abc = make_simplex({1, 2, 3});
  vector<Simplex> subsxs;
  abc.each_subsimplex(1, [&](const Simplex& s) { subsxs.push_back(s); });
  DegenerateSimplices dsxs(subsxs.begin(), subsxs.end());
  BOOST_CHECK(dsxs.is_collapsed(subsxs[0]));
  BOOST_CHECK(!dsxs.is_collapsed(abc));
}

BOOST_AUTO_TEST_CASE(IsAcyclicTests) {
  Simplex abc = make_simplex({1, 2, 3});
  vector<Simplex> subsxs;
  abc.each_subsimplex(1, [&](const Simplex& s) { subsxs.push_back(s); });
  DegenerateSimplices dsxs(subsxs.begin(), subsxs.end());
  BOOST_CHECK(!dsxs.is_acyclic(abc));
}

BOOST_AUTO_TEST_SUITE_END();

