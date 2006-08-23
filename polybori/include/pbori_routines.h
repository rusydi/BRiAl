// -*- c++ -*-
//*****************************************************************************
/** @file pbori_routines.h
 *
 * @author Alexander Dreyer
 * @date 2006-08-23
 *
 * This file includes files, which define function templates for internal use 
 * in the polybori library.
 *
 * @par Copyright:
 *   (c) 2006 by
 *   Dep. of Mathematics, Kaiserslautern University of Technology and @n
 *   Fraunhofer Institute for Industrial Mathematics (ITWM)
 *   D-67663 Kaiserslautern, Germany
 *
 * @internal 
 * @version \$Id$
 *
 * @par History:
 * @verbatim
 * $Log$
 * Revision 1.1  2006/08/23 14:24:54  dreyer
 * ADD: BooleSet::usedVariables and infrastructure
 *
 * @endverbatim
**/
//*****************************************************************************

// include basic definitions
#include "pbori_defs.h"

// include polybori algorithms and functionals
#include "pbori_algo.h"
#include "pbori_func.h"

#ifndef PBORI_ROUTINES_H_
#define PBORI_ROUTINES_H_


// Get routines, which add features related to decision diagrams
#include "pbori_routines_dd.h"

// Get routines, which add features related to Cudd library
#include "pbori_routines_cuddext.h"

#endif
