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

std::vector<BoolePolynomial> someNextDegreeSpolys(GroebnerStrategy& strat, std::size_t n)
{
    std::vector<BoolePolynomial> ret;

    strat.pairs.cleanTopByChainCriterion();
    auto sugar_degree = strat.pairs.queue.top().sugar;

    while(
            !strat.pairs.pairSetEmpty() &&
            strat.pairs.queue.top().sugar <= sugar_degree &&
            ret.size() < n
         )
    {
        ret.emplace_back(strat.nextSpoly());
        strat.pairs.cleanTopByChainCriterion();
    }

    return ret;
}

std::vector<BoolePolynomial> easy_linear_polynomials_func(const BoolePolynomial & p)
{
    std::vector<BoolePolynomial> ret;
    return ret;
}

void high_probability_polynomials_trick(const BoolePolynomial & p, GroebnerStrategy & strat)
{
    return;
}

void add_to_basis(GroebnerStrategy & strat, const BoolePolynomial & p,
        bool prot, bool full_prot, bool use_easy_linear_polynomials)
{
    if (p.isZero())
    {
        if (prot)
            std::cout << "-" << std::endl;
    }
    else
    {
        if (prot)
        {
            if (full_prot)
                std::cout << p << std::endl;
            std::cout << "Result: deg: " << p.deg() << "lm : " << p.lead() << "el: " << p.eliminationLength() << std::endl;
        }

        if (use_easy_linear_polynomials && p.leadDeg() > 2)
        {
            std::vector<BoolePolynomial> lin = easy_linear_polynomials_func(p);
            for(const auto & q : lin)
            {
                strat.addGeneratorDelayed(q);
            }
        }

        std::size_t old_len = strat.generators.size();
        strat.addAsYouWish(p);
        std::size_t new_len = strat.generators.size();

        if (new_len == 1 + old_len)
            high_probability_polynomials_trick(p, strat);

        if (prot)
            std::cout << "Generators : " << strat.generators.size() << std::endl;
    }
}

std::vector<BoolePolynomial> symmGB_F2_python(
    const std::vector<BoolePolynomial> & F,
    uint64_t    deg_bound,
    std::size_t over_deg_bound,
    bool        use_faugere,
    bool        use_noro,
    bool        opt_lazy,
    bool        opt_red_tail,
    double      max_growth,
    double      step_factor,
    bool        implications,
    bool        prot,
    bool        full_prot,
    std::size_t selection_size,
    bool        opt_exchange,
    bool        opt_allow_recursion,
    bool        ll,
    bool        opt_linear_algebra_in_last_block,
    std::size_t max_generators,
    bool        red_tail_deg_growth,
    std::string matrix_prefix,
    bool        modified_linear_algebra,
    bool        draw_matrices,
    bool        easy_linear_polynomials
)
{
    PBORI_ASSERT(!F.empty());

    if (use_noro && use_faugere)
    {
        throw std::runtime_error("symmGB_F2_python() : both use_noro and use_faugere are specified");
    }

    std::vector<BoolePolynomial> ret;

    GroebnerStrategy GBStrategy(F[0].ring());
    GBStrategy.generators.optRedTail = opt_red_tail;
    GBStrategy.optLazy = opt_lazy;
    GBStrategy.optExchange = opt_exchange;
    GBStrategy.optAllowRecursion = opt_allow_recursion;
    GBStrategy.enabledLog = prot;
    GBStrategy.generators.optLL = ll;
    GBStrategy.optModifiedLinearAlgebra = modified_linear_algebra;
    GBStrategy.optLinearAlgebraInLastBlock = opt_linear_algebra_in_last_block;
    GBStrategy.reduceByTailReduced = false;
    GBStrategy.generators.optRedTailDegGrowth = red_tail_deg_growth;

    GBStrategy.optDrawMatrices = draw_matrices;
    GBStrategy.matrixPrefix = matrix_prefix;

    for(const auto & f : F)
    {
        if (!f.isZero())
            GBStrategy.addGeneratorDelayed(f);
    }

    if (prot)
        std::cout << "added delayed" << std::endl;

    std::size_t loopcounter = 0;

    while(GBStrategy.pairs.queue.size() > 0)
    {
        if (max_generators && GBStrategy.generators.size() > max_generators)
        {
            throw std::length_error("symmGB_F2_python() : the number of generators exceeds max_generators");
        }
        loopcounter += 1;

        auto sugar_degree = GBStrategy.pairs.queue.top().sugar;

        if (prot)
        {
            std::cout << "Current loop : " << loopcounter << std::endl;
            std::cout << "Current Degree : " << sugar_degree << std::endl;
            std::cout << "Current Generator size : " << GBStrategy.generators.size() << std::endl;
        }

        if (sugar_degree > deg_bound && over_deg_bound == 0)
        {
            return GBStrategy.minimalizeAndTailReduce();
        }

        std::vector<BoolePolynomial> ps;
        if (sugar_degree > deg_bound)
        {
            ps = someNextDegreeSpolys(GBStrategy, over_deg_bound);
            over_deg_bound -= ps.size();
        }
        else
        {
            ps = someNextDegreeSpolys(GBStrategy, selection_size);
        }

        if (!ps.empty() && ps[0].ring().ordering().isDegreeOrder())
        {
            for(auto & p : ps)
                p = cheap_reductions(GBStrategy.generators, p);

            auto end = std::remove_if(ps.begin(), ps.end(),
                    [](const BoolePolynomial & p){ return p.isZero(); });
            ps.erase(end, ps.end());
            std::size_t min_degree;
            if (ps.size() > 0)
            {
                auto it = std::min_element(ps.begin(), ps.end(),
                        [](const BoolePolynomial & l, const BoolePolynomial & r)
                        {
                            return l.deg() < r.deg();
                        });
                std::size_t min_degree = it->deg();
            }

            std::vector<BoolePolynomial> new_ps;
            for(const auto & p : ps)
            {
                if (p.deg() <= min_degree)
                    new_ps.emplace_back(std::move(p));
                else
                    GBStrategy.addGeneratorDelayed(p);
            }
            ps = std::move(new_ps);
        }

        if (prot)
        {
            std::cout << "(" << GBStrategy.pairs.queue.size() << ")" << std::endl;
            std::cout << "start reducing" << std::endl;
            std::cout << "Chain Crit. : " << GBStrategy.chainCriterions << std::endl;
            std::cout << "VC: " << GBStrategy.variableChainCriterions << std::endl;
            std::cout << "EASYP " << GBStrategy.easyProductCriterions << std::endl;
            std::cout << "EXTP " << GBStrategy.extendedProductCriterions << std::endl;
            std::cout << ps.size() << " spolys added" << std::endl;
        }


        if (use_noro || use_faugere)
        {
            std::vector<BoolePolynomial> v;

            for(const auto & p : ps)
            {
                if (!p.isZero())
                {
                    v.push_back(p);
                }
            }

            if (use_noro)
                ret = GBStrategy.noroStep(v);
            else
                ret = GBStrategy.faugereStepDense(v);
        }
        else
        {
            std::vector<BoolePolynomial> v;
            for(const auto & p : ps)
            {
                if (!p.isZero())
                {
                    v.push_back(p);
                }
            }

            if (v.size() > 100)
            {
                ret = parallel_reduce(v, GBStrategy, int(step_factor*10), max_growth);
            }
            else
            {
                if (v.size() > 10)
                {
                    ret = parallel_reduce(v, GBStrategy, int(step_factor*30), max_growth);
                }
                else
                {
                    ret = parallel_reduce(v, GBStrategy, int(step_factor * 100), max_growth);
                }
            }

        }

        if (prot)
            std::cout << "end reducing" << std::endl;

        if (ret.size() > 0 && ret[0].ring().ordering().isDegreeOrder())
        {
            auto it = std::min_element(ret.begin(), ret.end(),
                    [](const BoolePolynomial & l, const BoolePolynomial & r)
                    {
                        return l.deg() < r.deg();
                    });
            std::size_t ret_min_deg = it->deg();
            std::vector<BoolePolynomial> new_ret;
            for(const auto & p : ret)
            {
                if (p.deg() == ret_min_deg)
                    new_ret.emplace_back(std::move(p));
                else
                    GBStrategy.addGeneratorDelayed(p);
            }
            ret = std::move(new_ret);
        }

        std::sort(ret.begin(), ret.end(),
                [](const BoolePolynomial & l, const BoolePolynomial & r)
                {
                    return l.lead() < r.lead();
                });

        std::vector<BoolePolynomial> res_cp(std::move(ret));

        for(const auto & p : res_cp)
        {
            std::size_t old_len = GBStrategy.generators.size();
            add_to_basis(GBStrategy, p, prot, full_prot, easy_linear_polynomials);
            if (implications && old_len == GBStrategy.generators.size() - 1)
            {
                GBStrategy.addNonTrivialImplicationsDelayed(
                        const_cast<const GroebnerStrategy&>(GBStrategy).generators[GBStrategy.generators.size() - 1]);
            }

            if (p.isOne())
            {
                if (prot)
                {
                    std::cout << "GB is 1" << std::endl;
                    return GBStrategy.minimalizeAndTailReduce();
                }
            }

            if (prot)
                std::cout << "(" << GBStrategy.pairs.queue.size() << ")" << std::endl;
        }
        GBStrategy.cleanTopByChainCriterion();
    }

    ret = GBStrategy.minimalizeAndTailReduce();
    return ret;
}



END_NAMESPACE_PBORIGB

