/*
 *
 */

#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>
using boost::test_tools::output_test_stream;

#include <polybori/groebner/symmGB.h>

USING_NAMESPACE_PBORI
USING_NAMESPACE_PBORIGB

BOOST_AUTO_TEST_SUITE(symmGBTestSuite)

BOOST_AUTO_TEST_CASE(TestSymmGB_F2_C)
{
    BoolePolyRing ring(5, BoolePolyRing::ordercodes::dp_asc);
    BooleVariable x0(4, ring), x1(3, ring), x2(2, ring), x3(1, ring), x4(0, ring);

    BoolePolynomial f0 = x3*x0 + x3*x1 + x4*x1 + x4*x3 + x0 + x1 + x2 + x3 + x4;
    BoolePolynomial f1 = x1*x0 + x2*x1 + x3*x1 + x3*x2 + x4*x2 + x4*x3 + x0;
    BoolePolynomial f2 = x2*x0 + x2*x1 + x3*x0 + x3*x1 + x4*x0 + x4*x1 + x1 + x2 + x3;
    BoolePolynomial f3 = x2*x1 + x3*x2 + x0 + x3 + 1;
    BoolePolynomial f4 = x2*x1 + x3*x1 + x4*x0 + x4*x2 + x4*x3 + x1 + x3 + 1;
    BoolePolynomial f5 = x2*x0 + x2*x1 + x3*x2 + x4*x2 + x3 + x4 + 1;
    BoolePolynomial f6 = x2*x0 + x3*x1 + x4*x1 + x0 + 1;

    std::vector<BoolePolynomial> F = {f0, f1, f2, f3, f4, f5, f6};

    BoolePolynomial g0 = x0;
    BoolePolynomial g1 = x1 + 1;
    BoolePolynomial g2 = x2 + 1;
    BoolePolynomial g3 = x3;
    BoolePolynomial g4 = x4 + 1;

    std::vector<BoolePolynomial> ans = {g0, g1, g2, g3, g4};
    std::sort(ans.begin(), ans.end(),
            [](const BoolePolynomial & f, const BoolePolynomial & g)
            {
                return f.lead() < g.lead();
            });

    std::vector<BoolePolynomial> G = symmGB_F2_C(F);

    BOOST_CHECK_EQUAL_COLLECTIONS(G.begin(), G.end(), ans.begin(), ans.end());
}

BOOST_AUTO_TEST_SUITE_END()
