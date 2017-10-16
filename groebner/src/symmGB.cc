/*
 * symmGB.cc
 *
 *  Created on: Oct 16, 2017
 *      Author: makarim
 */

#include <polybori/groebner/symmGB.h>

BEGIN_NAMESPACE_PBORIGB

std::vector<BoolePolynomial> symmGB_F2_C(
    const std::vector<BoolePolynomial> & F,
    bool        opt_exchange,
    uint64_t    deg_bound,
    bool        opt_lazy,
    std::size_t over_deg_bound,
    bool        opt_red_tail,
    double      max_growth,
    double      step_factor,
    bool        implications,
    bool        prot,
    bool        full_prot,
    std::size_t selection_size,
    bool        opt_allow_recursion,
    bool        use_noro,
    bool        use_faugere,
    bool        ll,
    bool        opt_linear_algebra_in_last_block,
    std::size_t max_generators,
    bool        red_tail_deg_growth,
    bool        modified_linear_algebra,
    std::string matrix_prefix,
    bool        draw_matrices
)
{
    PBORI_ASSERT(!F.empty());

    if (use_noro)
    {
        throw std::logic_error("noro not implemented for symmgb");
    }

    GroebnerStrategy GBStrategy(F[0].ring());
    GBStrategy.generators.optRedTail = opt_red_tail;
    GBStrategy.enabledLog = prot;
    GBStrategy.optLazy = opt_lazy;
    GBStrategy.optExchange = opt_exchange;
    GBStrategy.generators.optLL = ll;
    GBStrategy.optAllowRecursion = opt_allow_recursion;
    GBStrategy.optLinearAlgebraInLastBlock = opt_linear_algebra_in_last_block;
    GBStrategy.optModifiedLinearAlgebra = modified_linear_algebra;
    GBStrategy.matrixPrefix = matrix_prefix;
    GBStrategy.optDrawMatrices = draw_matrices;
    GBStrategy.generators.optRedTailDegGrowth = red_tail_deg_growth;

    GBStrategy.reduceByTailReduced = false;

    for(const auto & f : F)
    {
        if (!f.isZero())
            GBStrategy.addGeneratorDelayed(f);
    }
    GBStrategy.symmGB_F2();

    std::vector<BoolePolynomial> ret = GBStrategy.minimalizeAndTailReduce();

    return ret;
}

END_NAMESPACE_PBORIGB

