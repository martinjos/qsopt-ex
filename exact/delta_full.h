#ifndef __DELTA_FULL_H__
#define __DELTA_FULL_H__

#include "basicdefs.h"
#include "qstruct_mpq.h"

/** @file
 *  @ingroup DeltaSolver */
/** @defgroup DeltaSolver delta-complete solver
 *  @{ */

/* ========================================================================= */
/** @brief Callback to receive partial results from delta-complete full LP
 * method. */
typedef void (*delta_full_callback_t)(const mpq_QSdata * const p_mpq,
                                      const mpq_t * const x,
                                      const mpq_t * const y,
                                      const mpq_t obj_lo,
                                      const mpq_t obj_up,
                                      const mpq_t diff,
                                      const mpq_t delta,
                                      void *data);

/* ========================================================================= */
/** @brief Given an mpq_QSdata problem, solve the corresponding
 * delta-optimality problem exactly.
 * @param p_mpq problem for which to determine delta-optimality exactly.
 * @param delta the delta to use for determining delta-optimality; the maximum
 * distance to optimality.
 * @param x we store here a primal feasible solution to the problem (if primal
 * feasibility established).
 * @param y we store here a dual feasible solution to the problem (if dual
 * feasibility established), or a certificate of infeasibility for the problem
 * (if infeasibility established).
 * @param obj_lo we store here a lower bound on the objective (if established).
 * @param obj_up we store here an upper bound on the objective (if established).
 * @param ebasis if not null, use the given basis to start the iteration of
 * simplex, and store here the final basis (where applicable).
 * @param simplexalgo whether to use primal or dual simplex while solving the
 * delta-optimality problem.
 * @param status pointer to the integer where we will return the status of the
 * problem, either optimal, delta-optimal, unbounded, or infeasible (we could
 * also return time out).
 * @param delta_callback if not null, will be called if a delta-satisfying
 * result is found for some value greater than delta.
 * @param callback_data additional parameter to be passed to delta_callback.
 * @return zero on success, non-zero otherwise. */
int QSdelta_full_solver (mpq_QSdata * p_mpq,
                         const mpq_t delta,
                         mpq_t * const x,
                         mpq_t * const y,
                         mpq_t obj_lo,
                         mpq_t obj_up,
                         QSbasis * const ebasis,
                         int simplexalgo,
                         int *status,
                         delta_full_callback_t delta_callback,
                         void *callback_data);

/** @} */

#endif
