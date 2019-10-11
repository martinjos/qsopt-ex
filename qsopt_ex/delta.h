#ifndef __DELTA_H__
#define __DELTA_H__

#include "qstruct_mpq.h"

/** @file
 *  @ingroup DeltaSolver */
/** @defgroup DeltaSolver delta-complete solver
 *  @{ */

/* ========================================================================= */
/** @brief Given an mpq_QSdata problem, solve the corresponding
 * delta-feasibility problem exactly.
 * @param p_mpq problem for which to determine delta-feasibility exactly.
 * @param delta the delta to use for determining delta-feasibility; the maximum
 * perturbation of RHS/bounds required to make a delta-feasible solution
 * feasible; updated with the actual infeasibility if delta-feasibility
 * determined.
 * @param x if not null, we store here a delta-feasible solution to the problem
 * (if delta-feasibility established).
 * @param y if not null, we store here a certificate of infeasibility for the
 * problem (if infeasibility established).
 * @param ebasis if not null, use the given basis to start the iteration of
 * simplex, and store here the final basis (where applicable).
 * @param simplexalgo whether to use primal or dual simplex while solving the
 * delta-feasibility problem.
 * @param status pointer to the integer where we will return the status of the
 * problem, either feasible, delta-feasible or infeasible (we could also return
 * time out).
 * @return zero on success, non-zero otherwise. */
int QSdelta_solver (mpq_QSdata * p_orig,
                    mpq_t delta,
                    mpq_t * const x,
                    mpq_t * const y,
                    QSbasis * const ebasis,
                    int simplexalgo,
                    int *status);

/** @} */

#endif
