/* ========================================================================= */
/* ESolver "Exact Mixed Integer Linear Solver" provides some basic structures
 * and algorithms commons in solving MIP's
 *
 * Copyright (C) 2005, 2019 Daniel Espinoza and Martin Sidaway.
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
#ifndef __EXACT_DELTA_H__
#define __EXACT_DELTA_H__

#include "exact.h"

/* ========================================================================= */
/** @defgroup Esolver Esolver
 * Here we define an interface to solve the delta-satisfiability problem for
 * LP's exactly (#QSexact_delta_solver).
 * */
/** @file
 * @ingroup Esolver */
/** @addtogroup Esolver */
/** @{ */
/* ========================================================================= */

/* ========================================================================= */
int QSexact_delta_optimal_test (mpq_QSdata * p,
																mpq_t * p_sol,
																mpq_t * d_sol,
																QSbasis * basis,
																mpq_t delta,
																delta_callback_t delta_callback,
																mpq_t last_infeas,
																void *callback_data);

/* ========================================================================= */
/** @brief Given an mpq_QSdata problem, solve the delta-satisfiability problem
 * exactly.  The problem must be a feasible and bounded feasibility LP with
 * artificial variables for all constraints.
 * @param x if not null, we store here the primal solution to the 
 * problem (if it exist).
 * @param y if not null, we store here the dual solution to the
 * problem, 
 * @param p_mpq problem to solve for delta-satisfiability exactly.
 * @param status pointer to the integer where we will return the status
 * of the problem, either optimal, infeasible, or unbounded (we could also 
 * return time out).
 * @param simplexalgo whether to use primal or dual simplex while solving
 * to optimality the problem.
 * @param basis if not null, use the given basis to start the
 * iteration of simplex, and store here the optimal basis (if found).
 * @param delta the elementwise maximum infeasibility to accept in a
 * delta-satisfying assignment.
 * @param delta_callback if not null, will be called if a delta-satisfying
 * result is found for some value greater than delta.
 * @param callback_data additional parameter to be passed to delta_callback.
 * @return zero on success, non-zero otherwise. */
int QSexact_delta_solver (mpq_QSdata * p_mpq,
													mpq_t * const x,
													mpq_t * const y,
													QSbasis * const basis,
													int simplexalgo,
													int *status,
													mpq_t delta,
													delta_callback_t delta_callback,
													void *callback_data);

/* ========================================================================= */
/** @brief create sum-of-infeasibilities objective function for QSexact_delta_solver */
int QSexact_delta_create_soi_obj (mpq_QSdata *p_mpq);

/** @} */
/* ========================================================================= */
/* end of exact_delta.h */
#endif

