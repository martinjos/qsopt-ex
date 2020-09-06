/* ========================================================================= */
/* ESolver "Exact Mixed Integer Linear Solver" provides some basic structures
 * and algorithms commons in solving MIP's
 *
 * Copyright (C) 2005 Daniel Espinoza.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 * */
/* ========================================================================= */
#ifndef __EXACT_BASIS_H__
#define __EXACT_BASIS_H__

#include <gmp.h>

#include "basicdefs.h"
#include "eg_io.h"
#include "eg_lpnum.h"
#include "eg_lpnum.dbl.h"
#include "eg_lpnum.mpq.h"
#include "eg_lpnum.mpf.h"
#include "lpdata_dbl.h"
#include "lpdata_mpq.h"
#include "lpdata_mpf.h"
#include "qsopt_mpf.h" /* mpf_QSset_precision */
#include "qstruct_dbl.h"
#include "qstruct_mpq.h"
#include "qstruct_mpf.h"

/* ========================================================================= */
/** @defgroup Esolver Esolver
 * */
/** @file
 * @ingroup Esolver */
/** @addtogroup Esolver */
/** @{ */
/* ========================================================================= */

/* ========================================================================= */
/** @brief get the status for a given basis in rational arithmetic, it should
 * also leave everything set to get primal/dual solutions when needed.
 * @param p_mpq       the problem data.
 * @param status      where to store status of basis (optimal, infeasible,
 *                    unbounded, or unsolved).
 * @param basis       basis to be tested.
 * @param msg_lvl     message level.
 * @param simplexalgo may be updated with PRIMAL_SIMPLEX to indicate that
 *                    subsequent solves should try primal.
 */
int QSexact_basis_status (mpq_QSdata * p_mpq,
                          int *status,
                          QSbasis * const basis,
                          const int msg_lvl,
                          int *const simplexalgo);

/* ========================================================================= */
/** @brief get the status for a given basis in rational arithmetic, it should
 * also leave everything set to get primal/dual solutions when needed.
 * @param p_mpq       the problem data.
 * @param status      where to store status of basis (optimal, infeasible,
 *                    unbounded, or unsolved).
 * @param basis       basis to be tested.
 * @param msg_lvl     message level.
 * @param simplexalgo may be updated with PRIMAL_SIMPLEX to indicate that
 *                    subsequent solves should try primal.
 * */
int QSdelta_basis_status (mpq_QSdata * p_mpq,
                          int *status,
                          QSbasis * const basis,
                          const int msg_lvl,
                          int *const simplexalgo);

/* ========================================================================= */
/** @brief get the status for a given basis in rational arithmetic, it should
 * also leave everything set to get primal/dual solutions when needed.
 * @param p_mpq       the problem data.
 * @param status      where to store status of basis (optimal, infeasible,
 *                    unbounded, or unsolved).
 * @param basis       basis to be tested.
 * @param msg_lvl     message level.
 * @param simplexalgo may be updated with PRIMAL_SIMPLEX to indicate that
 *                    subsequent solves should try primal.
 * */
int QSdelta_full_basis_status (mpq_QSdata * p_mpq,
                               int *status,
                               QSbasis * const basis,
                               const int msg_lvl,
                               int *const simplexalgo);

/** @} */
/* ========================================================================= */
/* end of exact_basis.h */
#endif

