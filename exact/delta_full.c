/** @file
 *  @ingroup DeltaSolver */
/** @addtogroup DeltaSolver
 *  @{ */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "delta_full.h"
#include "delta.h"
#include "exact_basis.h"

#include "qstruct_dbl.h"
#include "qstruct_mpf.h"

#include "basis_mpq.h"
#include "exact.h"
#include "except.h"
#include "editor_dbl.h"
#include "editor_mpf.h"
#include "eg_macros.h"
#include "eg_timer.h"
#include "eg_numutil_mpq.h"
#include "fct_mpq.h"
#include "qsopt_mpf.h"
#include "qsopt_dbl.h"
#include "qsopt_mpq.h"
#include "simplex_mpq.h"
#include "dump.h"

#include <assert.h>

// Source size is known, because it was allocated using __EGlpNumAllocArray.
// Target size is not known - this must be part of the function contract.
void copy_array (mpq_t * const target,
                 const mpq_t * const source)
{
  unsigned sz = __EGlpNumArraySize (source);
  while (sz--)
    mpq_set (target[sz], source[sz]);
}

/* ========================================================================= */
/** @brief print into screen (if enable) a message indicating that we have
 * successfully prove infeasibility, and save (if y is non
 * NULL ) the dual ray solution provided in y_mpq.
 * @param p_mpq the problem data.
 * @param y where to store the optimal dual solution (if not null).
 * @param y_mpq  the optimal dual solution.
 * */
/* ========================================================================= */
static void infeasible_output (mpq_QSdata * p_mpq,
                               mpq_t * const y,
                               const mpq_t * const y_mpq)
{
  if (p_mpq->simplex_display)
  {
    QSlog("Problem is infeasible");
  }
  if (y)
  {
    copy_array (y, y_mpq);
  }
}

/* ========================================================================= */
/** @brief print into screen (if enable) a message indicating that we have
 * successfully proven unboundedness.
 * @param p_mpq the problem data.
 * */
/* ========================================================================= */
static void unbounded_output (mpq_QSdata * p_mpq)
{
  if (p_mpq->simplex_display)
  {
    QSlog("Problem is unbounded");
  }
}

/* ========================================================================= */
/** @brief print into screen (if enable) a message indicating that we have
 * successfully proven optimality.
 * @param p_mpq the problem data.
 * */
/* ========================================================================= */
static void optimal_output (mpq_QSdata * p_mpq)
{
  if (p_mpq->simplex_display)
  {
    QSlog("Found optimal solution");
  }
}

/* ========================================================================= */
/** @brief print into screen (if enable) a message indicating that we have
 * successfully proven delta-optimality.
 * @param p_mpq the problem data.
 * */
/* ========================================================================= */
static void delta_optimal_output (mpq_QSdata * p_mpq)
{
  if (p_mpq->simplex_display)
  {
    QSlog("Found delta-optimal solution");
  }
}

/* ========================================================================= */
/** @brief copy the dual solution out of p_mpq.
 * @param y where to store the feasible dual solution.
 * @param p_mpq the problem data.
 * */
/* ========================================================================= */
static int copy_y (mpq_t * const y,
                   const mpq_QSdata * p_mpq)
{
  int rval = 0;
  int i, col;
  mpq_lpinfo* lp = p_mpq->lp;

  // QSexact_basis_status() has already been called, so we have the basis,
  // and the y value (piz).
  for (i = 0; i < lp->nrows; i++)
  {
    mpq_EGlpNumCopy (y[i], lp->piz[i]);
  }

CLEANUP:

  EG_RETURN (rval);
}

// There is simply no way to use mpq___EGlpNumInnProd, because it requires
// result to be an mpq_t* (__mpq_struct(*)[1]), which is impossible starting
// from an mpq_t (__mpq_struct[1]) parameter (__mpq_struct*).
static void my_inner_prod (mpq_t result, mpq_t *const a, mpq_t *const b, const size_t len)
{
  mpq_EGlpNumZero(result);
  for (size_t i = len; i--;)
  {
    mpq_EGlpNumAddInnProdTo(result, a[i], b[i]);
  }
}

static int judge_basis (int *judgement,
                        int retry_get_infeas_array,
                        mpq_QSdata * p_mpq,
                        int *status,
                        QSbasis * const basis,
                        const int msg_lvl,
                        const mpq_t delta,
                        mpq_t * const x,
                        mpq_t * const y,
                        mpq_t obj_lo,
                        mpq_t obj_up,
                        mpq_t *obj_coefs,
                        mpq_t *rhs_coefs,
                        int *have_primal,
                        int *have_dual)
{
  int rval = 0;
  mpq_t *y_mpq = 0;
  mpq_t diff;

  mpq_EGlpNumInitVar (diff);
  *judgement = 0;

  MESSAGE (msg_lvl, "Basis hash is 0x%016lX", QSexact_basis_hash(basis));
  EGcallD(QSexact_basis_status (p_mpq, status, basis, msg_lvl, NULL));
  if (QS_LP_INFEASIBLE == *status)
  {
    y_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->nrows);
    int mpq_QSget_infeas_array_rval = mpq_QSget_infeas_array (p_mpq, y_mpq);
    if (mpq_QSget_infeas_array_rval)
    {
      if (retry_get_infeas_array)
      {
        *judgement = 2;
        goto CLEANUP;
      }
      else
        EGcallD (mpq_QSget_infeas_array_rval);
    }
    infeasible_output (p_mpq, y, y_mpq);
    *judgement = 1;
    goto CLEANUP;
  }
  else if (QS_LP_UNBOUNDED == *status)
  {
    // TODO: return certificate of unboundedness?
    unbounded_output (p_mpq);
    *judgement = 1;
    goto CLEANUP;
  }
  else
  {
    if (QS_LP_OPTIMAL == *status || p_mpq->lp->probstat.primal_feasible)
    {
      assert (p_mpq->lp->basisstat.primal_feasible);
      EGcallD (QSdelta_copy_x (x, p_mpq));
      mpq_t primal_obj;
      mpq_init (primal_obj);
      my_inner_prod (primal_obj, obj_coefs, x, p_mpq->qslp->nstruct);
#ifndef NDEBUG
      mpq_ILLfct_compute_pobj (p_mpq->lp);
      if (mpq_cmp (p_mpq->lp->pobjval, primal_obj) != 0)
      {
        MESSAGE (msg_lvl, "Oops: p_mpq->lp->pobjval = %lf but primal_obj = %lf",
                 mpq_get_d (p_mpq->lp->pobjval),
                 mpq_get_d (primal_obj));
      }
      assert (mpq_cmp (p_mpq->lp->pobjval, primal_obj) == 0);
#endif
      if (!*have_primal ||
          (p_mpq->qslp->objsense == QS_MIN ? mpq_cmp (primal_obj, obj_up) < 0
                                           : mpq_cmp (primal_obj, obj_lo) > 0))
      {
        mpq_set (p_mpq->qslp->objsense == QS_MIN ? obj_up : obj_lo,
                 primal_obj);
        MESSAGE (msg_lvl, "Primal feasible: set %s to %lf",
                 p_mpq->qslp->objsense == QS_MIN ? "obj_up" : "obj_lo",
                 mpq_get_d (primal_obj));
      }
      mpq_clear (primal_obj);
      *have_primal = 1;
    }
    if (QS_LP_OPTIMAL == *status || p_mpq->lp->probstat.dual_feasible)
    {
      assert (p_mpq->lp->basisstat.dual_feasible);
      EGcallD (copy_y (y, p_mpq));
      mpq_ILLfct_compute_dobj (p_mpq->lp);
      if (!*have_dual ||
          (p_mpq->qslp->objsense == QS_MIN
            ? mpq_cmp (p_mpq->lp->dobjval, obj_lo) > 0
            : mpq_cmp (p_mpq->lp->dobjval, obj_up) < 0))
      {
        mpq_set (p_mpq->qslp->objsense == QS_MIN ? obj_lo : obj_up,
                 p_mpq->lp->dobjval);
        MESSAGE (msg_lvl, "Dual feasible: set %s to %lf",
                 p_mpq->qslp->objsense == QS_MIN ? "obj_lo" : "obj_up",
                 mpq_get_d (p_mpq->lp->dobjval));
        QSlog ("bz:");
        mpq_QSdump_bz (p_mpq);
        QSlog ("piz:");
        mpq_QSdump_piz (p_mpq);
      }
      *have_dual = 1;
    }
#if 0
    if (QS_LP_OPTIMAL == *status)
    {
      // Shortcut: can avoid difference check
      optimal_output (p_mpq);
      *judgement = 1;
      goto CLEANUP;
    }
#endif
    if (*have_primal && *have_dual)
    {
      // Could be optimal or delta-optimal
      mpq_sub (diff, obj_up, obj_lo);
      assert (mpq_sgn (diff) >= 0);
      if (mpq_sgn (diff) == 0)
      {
        // Optimal
        optimal_output (p_mpq);
        *status = QS_LP_OPTIMAL;
        *judgement = 1;
        goto CLEANUP;
      }
      else if (mpq_cmp (diff, delta) <= 0)
      {
        // Delta-optimal
        delta_optimal_output (p_mpq);
        *status = QS_LP_DELTA_OPTIMAL;
        *judgement = 1;
        goto CLEANUP;
      }
    }
  }

CLEANUP:
  mpq_EGlpNumFreeArray (y_mpq);
  mpq_EGlpNumClearVar (diff);
  EG_RETURN (rval);  
}

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
 * @return zero on success, non-zero otherwise. */
int QSdelta_full_solver (mpq_QSdata * p_mpq,
                         const mpq_t delta,
                         mpq_t * const x,
                         mpq_t * const y,
                         mpq_t obj_lo,
                         mpq_t obj_up,
                         QSbasis * const ebasis,
                         int simplexalgo,
                         int *status)
{
  /* local variables */
  int last_status = 0, last_iter = 0;
  QSbasis *basis = 0;
  unsigned precision = EGLPNUM_PRECISION;
  int rval = 0,
    it = QS_EXACT_MAX_ITER;
  dbl_QSdata *p_dbl = 0;
  mpf_QSdata *p_mpf = 0;
  double *x_dbl = 0,
   *y_dbl = 0;
  mpf_t *x_mpf = 0,
   *y_mpf = 0;
  int have_primal = 0, have_dual = 0;
  mpq_t *obj_coefs = 0,
        *rhs_coefs = 0;
  mpq_t zero;

  *status = QS_LP_UNSOLVED;
  mpq_EGlpNumInitVar (zero);  // Inits to zero

  if (!p_mpq || !p_mpq->qslp)
  {
    QSlog("LP data not available");
    rval = 1;
    goto CLEANUP;
  }
  int const msg_lvl = __QS_SB_VERB <= DEBUG ? 0: (1 - p_mpq->simplex_display) * 10000;
  mpq_set_si(obj_lo, 0, 1);  // Result is exact unless returning actual
  mpq_set_si(obj_up, 0, 1);  // objective function value
  /* save the problem if we are really debugging */
  if(DEBUG >= __QS_SB_VERB)
  {
    EGcallD(mpq_QSwrite_prob(p_mpq, "qsxprob.lp","LP"));
  }
  /* Get objective function and RHS */
  obj_coefs = mpq_EGlpNumAllocArray (p_mpq->qslp->nstruct);
  EGcallD (mpq_QSget_obj (p_mpq, obj_coefs));
  rhs_coefs = mpq_EGlpNumAllocArray (p_mpq->qslp->nrows);
  EGcallD (mpq_QSget_rhs (p_mpq, rhs_coefs));
  /* try first with doubles */
  if (p_mpq->simplex_display || DEBUG >= __QS_SB_VERB)
  {
    QSlog("Trying double precision");
  }
  p_dbl = QScopy_prob_mpq_dbl (p_mpq, "dbl_problem");
  if(__QS_SB_VERB <= DEBUG && !p_dbl->simplex_display) p_dbl->simplex_display = 1;
  if (ebasis && ebasis->nstruct)
    dbl_QSload_basis (p_dbl, ebasis);
  if (dbl_ILLeditor_solve (p_dbl, simplexalgo))
  {
    MESSAGE(p_mpq->simplex_display ? 0: __QS_SB_VERB,
            "double approximation failed, code %d, "
            "continuing in extended precision", rval);
    goto MPF_PRECISION;
  }
  EGcallD(dbl_QSget_status (p_dbl, status));
  if ((*status == QS_LP_INFEASIBLE) &&
      (p_dbl->lp->final_phase != PRIMAL_PHASEI) &&
      (p_dbl->lp->final_phase != DUAL_PHASEII))
    dbl_QSopt_primal (p_dbl, status);
  EGcallD(dbl_QSget_status (p_dbl, status));
  last_status = *status;
  EGcallD(dbl_QSget_itcnt(p_dbl, 0, 0, 0, 0, &last_iter));
  /* deal with the problem depending on what status we got from our optimizer */
  if (QS_LP_OPTIMAL == *status || QS_LP_UNBOUNDED == *status || QS_LP_INFEASIBLE == *status)
  {
    basis = dbl_QSget_basis (p_dbl);
    int judgement = 0;
    EGcallD (judge_basis (&judgement, 1, p_mpq, status, basis, msg_lvl,
                          delta, x, y, obj_lo, obj_up, obj_coefs, rhs_coefs,
                          &have_primal, &have_dual));
    MESSAGE(msg_lvl, "judge_basis returned with judgement = %d, *status = %d",
            judgement, *status);
    if (judgement == 2)
    {
      MESSAGE(p_mpq->simplex_display ? 0 : __QS_SB_VERB, "double approximation"
              " failed, code %d, continuing in extended precision", rval);
      goto MPF_PRECISION;
    }
    else if (judgement)
    {
      MESSAGE(msg_lvl, "Basis judgement made, quitting");
      goto CLEANUP;
    }
  }
  else
  {
    MESSAGE (msg_lvl, "Floating-point solver reports failure; trying next precision");
  }
  IFMESSAGE(p_mpq->simplex_display,"Retrying in extended precision");
  /* if we reach this point, then we have to keep going, we use the previous
   * basis ONLY if the previous precision thinks that it has the optimal
   * solution, otherwise we start from scratch. */
  precision = 128;
  MPF_PRECISION:
  dbl_QSfree_prob (p_dbl);
  p_dbl = 0;
  /* try with multiple precision floating points */
  for (; it--; precision = (unsigned) (precision * 1.5))
  {
    QSexact_set_precision (precision);
    if (p_mpq->simplex_display || DEBUG >= __QS_SB_VERB)
    {
      QSlog("Trying mpf with %u bits", precision);
    }
    p_mpf = QScopy_prob_mpq_mpf (p_mpq, "mpf_problem");
    if(DEBUG >= __QS_SB_VERB)
    {
      EGcallD(mpf_QSwrite_prob(p_mpf, "qsxprob.mpf.lp","LP"));
    }
    if(__QS_SB_VERB <= DEBUG && !p_mpf->simplex_display) p_mpf->simplex_display = 1;
    simplexalgo = PRIMAL_SIMPLEX;
    if(!last_iter) last_status = QS_LP_UNSOLVED;
    if(last_status == QS_LP_OPTIMAL || last_status == QS_LP_INFEASIBLE)
    {
      if (p_mpq->simplex_display || DEBUG >= __QS_SB_VERB)
      {
        QSlog("Reusing previous basis");
      }
      if (basis)
      {
        EGcallD(mpf_QSload_basis (p_mpf, basis));
        mpf_QSfree_basis (basis);
        simplexalgo = DUAL_SIMPLEX;
        basis = 0;
      }
      else if (ebasis && ebasis->nstruct)
      {
        mpf_QSload_basis (p_mpf, ebasis);
        simplexalgo = DUAL_SIMPLEX;
      }
    }
    else
    {
      if(p_mpf->basis)
      {
        mpf_ILLlp_basis_free(p_mpf->basis);
        p_mpf->lp->basisid = -1;
        p_mpf->factorok = 0;
      }
      if (p_mpq->simplex_display || DEBUG >= __QS_SB_VERB)
      {
        QSlog("Not using previous basis");
      }
    }
    if (mpf_ILLeditor_solve (p_mpf, simplexalgo))
    {
      if (p_mpq->simplex_display || DEBUG >= __QS_SB_VERB)
      {
        QSlog("mpf_%u precision failed, error code %d, continuing with "
                    "next precision", precision, rval);
       }
      goto NEXT_PRECISION;
    }
    EGcallD(mpf_QSget_status (p_mpf, status));
    if ((*status == QS_LP_INFEASIBLE) &&
        (p_mpf->lp->final_phase != PRIMAL_PHASEI) &&
        (p_mpf->lp->final_phase != DUAL_PHASEII))
      mpf_QSopt_primal (p_mpf, status);
    EGcallD(mpf_QSget_status (p_mpf, status));
    last_status = *status;
    EGcallD(mpf_QSget_itcnt(p_mpf, 0, 0, 0, 0, &last_iter));
    /* deal with the problem depending on status we got from our optimizer */
    if (QS_LP_OPTIMAL == *status || QS_LP_UNBOUNDED == *status || QS_LP_INFEASIBLE == *status)
    {
      basis = mpf_QSget_basis (p_mpf);
      int judgement = 0;
      EGcallD (judge_basis (&judgement, 0, p_mpq, status, basis, msg_lvl,
                            delta, x, y, obj_lo, obj_up, obj_coefs, rhs_coefs,
                            &have_primal, &have_dual));
      MESSAGE(msg_lvl, "judge_basis returned with judgement = %d, *status = %d",
              judgement, *status);
      if (judgement) {
        MESSAGE(msg_lvl, "Basis judgement made, quitting");
        goto CLEANUP;
      }
    }
    else
    {
      MESSAGE (msg_lvl, "Floating-point solver reports failure; trying next precision");
    }
  NEXT_PRECISION:
    mpf_QSfree_prob (p_mpf);
    p_mpf = 0;
  }

  MESSAGE(msg_lvl, "Iteration limit reached");
  *status = QS_LP_ITER_LIMIT;

CLEANUP:
  dbl_EGlpNumFreeArray (x_dbl);
  dbl_EGlpNumFreeArray (y_dbl);
  mpq_EGlpNumFreeArray (obj_coefs);
  mpq_EGlpNumFreeArray (rhs_coefs);
  mpf_EGlpNumFreeArray (x_mpf);
  mpf_EGlpNumFreeArray (y_mpf);
  mpq_EGlpNumClearVar (zero);
  if (ebasis && basis)
  {
    ILL_IFFREE (ebasis->cstat, char);
    ILL_IFFREE (ebasis->rstat, char);
    ebasis->nstruct = basis->nstruct;
    ebasis->nrows = basis->nrows;
    ebasis->cstat = basis->cstat;
    ebasis->rstat = basis->rstat;
    basis->cstat = basis->rstat = 0;
  }
  mpq_QSfree_basis (basis);
  dbl_QSfree_prob (p_dbl);
  mpf_QSfree_prob (p_mpf);

  EG_RETURN (rval);
}


/** @} */
