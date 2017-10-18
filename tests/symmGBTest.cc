/*
 *
 */

#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>
using boost::test_tools::output_test_stream;

#include <polybori/groebner/symmGB.h>

USING_NAMESPACE_PBORI
USING_NAMESPACE_PBORIGB

BoolePolyRing R(5, BoolePolyRing::ordercodes::dp_asc);

struct I
{
    I() : ring(R), x0(4, R), x1(3, R), x2(2, R), x3(1, R), x4(0, R)
    {
        BoolePolynomial f0 = x3*x0 + x3*x1 + x4*x1 + x4*x3 + x0 + x1 + x2 + x3 + x4;
        BoolePolynomial f1 = x1*x0 + x2*x1 + x3*x1 + x3*x2 + x4*x2 + x4*x3 + x0;
        BoolePolynomial f2 = x2*x0 + x2*x1 + x3*x0 + x3*x1 + x4*x0 + x4*x1 + x1 + x2 + x3;
        BoolePolynomial f3 = x2*x1 + x3*x2 + x0 + x3 + 1;
        BoolePolynomial f4 = x2*x1 + x3*x1 + x4*x0 + x4*x2 + x4*x3 + x1 + x3 + 1;
        BoolePolynomial f5 = x2*x0 + x2*x1 + x3*x2 + x4*x2 + x3 + x4 + 1;
        BoolePolynomial f6 = x2*x0 + x3*x1 + x4*x1 + x0 + 1;

        generators.reserve(7);
        generators.emplace_back(std::move(f0));
        generators.emplace_back(std::move(f1));
        generators.emplace_back(std::move(f2));
        generators.emplace_back(std::move(f3));
        generators.emplace_back(std::move(f4));
        generators.emplace_back(std::move(f5));
        generators.emplace_back(std::move(f6));

        BoolePolynomial g0 = x0;
        BoolePolynomial g1 = x1 + 1;
        BoolePolynomial g2 = x2 + 1;
        BoolePolynomial g3 = x3;
        BoolePolynomial g4 = x4 + 1;

        gb.reserve(5);
        gb.emplace_back(std::move(g0));
        gb.emplace_back(std::move(g1));
        gb.emplace_back(std::move(g2));
        gb.emplace_back(std::move(g3));
        gb.emplace_back(std::move(g4));

        std::sort(gb.begin(), gb.end(),
            [](const BoolePolynomial & f, const BoolePolynomial & g)
            {
                return f.lead() < g.lead();
            });
    }

    BoolePolyRing ring;
    BooleVariable x0, x1, x2, x3, x4;
    std::vector<BoolePolynomial> generators, gb;
};

BOOST_FIXTURE_TEST_SUITE(TestSymmGB_F2, I)

BOOST_AUTO_TEST_CASE(TestSymmGB_F2_C)
{
    std::vector<BoolePolynomial> G = symmGB_F2_C(generators);
    BOOST_CHECK_EQUAL_COLLECTIONS(G.begin(), G.end(), gb.begin(), gb.end());
}

/*
BOOST_AUTO_TEST_CASE(TestSymmGB_F2_python)
{
    std::vector<BoolePolynomial> G = symmGB_F2_python(generators);
    BOOST_CHECK_EQUAL_COLLECTIONS(G.begin(), G.end(), gb.begin(), gb.end());
}
*/

BOOST_AUTO_TEST_SUITE_END()
