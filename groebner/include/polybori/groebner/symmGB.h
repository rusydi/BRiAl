#ifndef polybori_groebner_symmGB_h_
#define polybori_groebner_symmGB_h_

#include "groebner_defs.h"
#include "GroebnerStrategy.h"
#include "nf.h"

#include <vector>

BEGIN_NAMESPACE_PBORIGB

std::vector<BoolePolynomial> symmGB_F2_C(
    const std::vector<BoolePolynomial> & F,
    bool               opt_exchange=true,
    uint64_t           deg_bound=1000000000000,
    bool               opt_lazy=false,
    std::size_t        over_deg_bound=0,
    bool               opt_red_tail=true,
    double             max_growth=2.0,
    double             step_factor=1.0,
    bool               implications=false,
    bool               prot=false,
    bool               full_prot=false,
    std::size_t        selection_size=1000,
    bool               opt_allow_recursion=false,
    bool               use_noro=false,
    bool               use_faugere=false,
    bool               ll=false,
    bool               opt_linear_algebra_in_last_block=true,
    std::size_t        max_generators=0,
    bool               red_tail_deg_growth=true,
    bool               modified_linear_algebra=true,
    std::string        matrix_prefix="",
    bool               draw_matrices=false
);

std::vector<BoolePolynomial> symmGB_F2_python(
    const std::vector<BoolePolynomial> & F,
    uint64_t           deg_bound=1000000000000,
    std::size_t        over_deg_bound=0,
    bool               use_faugere=false,
    bool               use_noro=false,
    bool               opt_lazy=true,
    bool               opt_red_tail=true,
    double             max_growth=2.0,
    double             step_factor=1.0,
    bool               implications=false,
    bool               prot=true,
    bool               full_prot=false,
    std::size_t        selection_size=1000,
    bool               opt_exchange=true,
    bool               opt_allow_recursion=false,
    bool               ll=false,
    bool               opt_linear_algebra_in_last_block=true,
    std::size_t        max_generators=0,
    bool               red_tail_deg_growth=true,
    std::string        matrix_prefix="mat",
    bool               modified_linear_algebra=true,
    bool               draw_matrices=false,
    bool               easy_linear_polynomials=false
);

std::vector<BoolePolynomial> someNextDegreeSpolys(
    GroebnerStrategy& strat,
    int n
);

std::vector<BoolePolynomial> easy_linear_polynomials_func(
    const BoolePolynomial & p
);

void high_probability_polynomials_trick(
    const BoolePolynomial & p,
    GroebnerStrategy & strat
);

END_NAMESPACE_PBORIGB

#endif
