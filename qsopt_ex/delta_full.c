/** @file
 *  @ingroup DeltaSolver */
/** @addtogroup DeltaSolver
 *  @{ */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "delta_full.h"

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

/* ========================================================================= */
/** @brief Used as separator while printing output to the screen (controlled by
 * enabling simplex_display in the mpq_QSdata */
/* ========================================================================= */
static const char __sp[81] =
  "================================================================================";

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
/** @brief copy the primal solution out of p_mpq.
 * @param x where to store the feasible primal solution.
 * @param p_mpq the problem data.
 * */
/* ========================================================================= */
static int copy_x (mpq_t * const x,
                   const mpq_QSdata * p_mpq)
{
  int rval = 0;
  int i, col;
  mpq_t *tempx = 0;
  mpq_lpinfo* lp = p_mpq->lp;
  mpq_ILLlpdata* qslp = p_mpq->qslp;

  // Populate x with values of structural variables
  if (lp->nrows != qslp->nrows ||
      lp->ncols != qslp->ncols ||
      lp->nnbasic != qslp->nstruct ||
      lp->ncols != lp->nrows + lp->nnbasic)
  {
    QSlog("Unexpected condition: lp and qslp dimensions do not match");
    rval = 1;
    ILL_CLEANUP;
  }
  tempx = mpq_EGlpNumAllocArray (lp->ncols);
  // Set basic variables
  for (i = 0; i < lp->nrows; i++)
    mpq_EGlpNumCopy (tempx[lp->baz[i]], lp->xbz[i]);
  // Set non-basic variables
  for (i = 0; i < lp->nnbasic; i++)
  {
    col = lp->nbaz[i];
    if (lp->vstat[col] == STAT_UPPER)
      mpq_EGlpNumCopy (tempx[col], lp->uz[col]);
    else if (lp->vstat[col] == STAT_LOWER)
      mpq_EGlpNumCopy (tempx[col], lp->lz[col]);
    else
      mpq_EGlpNumZero (tempx[col]);
  }
  // Get structural variables
  for (i = 0; i < qslp->nstruct; i++)
  {
    mpq_EGlpNumCopy (x[i], tempx[qslp->structmap[i]]);
  }

CLEANUP:

  mpq_EGlpNumFreeArray (tempx);

  EG_RETURN (rval);
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

  // QSdelta_full_basis_status() has already been called, so we have the basis,
  // which is all that is needed.
  lp->piz = mpq_EGlpNumAllocArray (lp->nrows);
  mpq_ILLfct_compute_piz (lp);
  for (i = 0; i < lp->nrows; i++)
  {
    mpq_EGlpNumCopy (y[i], lp->piz[i]);
  }

CLEANUP:

  EG_RETURN (rval);
}

/* ========================================================================= */
/** @brief print into screen (if enable) a message indicating that we have
 * successfully prove feasibility.
 * @param p_mpq the problem data.
 * */
/* ========================================================================= */
static int feasible_output (mpq_QSdata * p_mpq,
                            mpq_t * const x)
{
  int rval = 0;

  if (p_mpq->simplex_display)
    QSlog("Problem is feasible");
  if (x)
    EGcallD(copy_x (x, p_mpq));

CLEANUP:

  EG_RETURN (rval);
}

/* ========================================================================= */
/** @brief get the status for a given basis in rational arithmetic, it should
 * also leave everything set to get primal/dual solutions when needed.
 * */
int QSdelta_full_basis_status (mpq_QSdata * p_mpq,
                               int *status,
                               QSbasis * const basis,
                               const int msg_lvl,
                               int *const simplexalgo)
{
  int rval = 0,
  singular;
  mpq_feas_info fi;
  EGtimer_t local_timer;
  mpq_t zero;
  mpq_EGlpNumInitVar (zero);  // Inits to zero
  mpq_EGlpNumInitVar (fi.totinfeas);
  EGtimerReset (&local_timer);
  EGtimerStart (&local_timer);
  EGcallD(mpq_QSload_basis (p_mpq, basis));
  if (p_mpq->cache)
  {
    mpq_ILLlp_cache_free (p_mpq->cache);
    mpq_clear (p_mpq->cache->val);
    ILL_IFFREE (p_mpq->cache, mpq_ILLlp_cache);
  }
  p_mpq->qstatus = QS_LP_MODIFIED;
  if(p_mpq->qslp->sinfo)
  {
    mpq_ILLlp_sinfo_free(p_mpq->qslp->sinfo);
    ILL_IFFREE(p_mpq->qslp->sinfo, mpq_ILLlp_sinfo);
  }
  if(p_mpq->qslp->rA)
  {
    mpq_ILLlp_rows_clear (p_mpq->qslp->rA);
    ILL_IFFREE (p_mpq->qslp->rA, mpq_ILLlp_rows);
  }
  mpq_free_internal_lpinfo (p_mpq->lp);
  mpq_init_internal_lpinfo (p_mpq->lp);
  EGcallD(mpq_build_internal_lpinfo (p_mpq->lp));
  mpq_ILLfct_set_variable_type (p_mpq->lp);
  EGcallD(mpq_ILLbasis_load (p_mpq->lp, p_mpq->basis));
  EGcallD(mpq_ILLbasis_factor (p_mpq->lp, &singular));
  memset (&(p_mpq->lp->probstat), 0, sizeof (mpq_lp_status_info));
  memset (&(p_mpq->lp->basisstat), 0, sizeof (mpq_lp_status_info));
  mpq_ILLfct_compute_xbz (p_mpq->lp);
  mpq_ILLfct_check_pfeasible (p_mpq->lp, &fi, mpq_zeroLpNum);
  p_mpq->lp->final_phase = PRIMAL_PHASEI;  // For mpq_QSget_infeas_array
  p_mpq->lp->pIpiz = mpq_EGlpNumAllocArray (p_mpq->lp->nrows);
  p_mpq->lp->pIdz = mpq_EGlpNumAllocArray (p_mpq->lp->nnbasic);
  mpq_ILLfct_compute_phaseI_piz (p_mpq->lp);
  mpq_ILLfct_compute_phaseI_dz (p_mpq->lp);

  if (p_mpq->simplex_display >= 2)
  {
    unsigned sz;
    if (p_mpq->simplex_display >= 3)
    {
      mpq_QSdump_prob(p_mpq);
      mpq_QSdump_basis(p_mpq);
    }
    QSlog("QSdelta_full_basis_status: xnbz =");
    mpq_QSdump_xnbz(p_mpq);
    QSlog("QSdelta_full_basis_status: xbz =");
    mpq_QSdump_xbz(p_mpq);
    QSlog("QSdelta_full_basis_status: bfeas =");
    mpq_QSdump_bfeas(p_mpq);
    QSlog("QSdelta_full_basis_status: pIpiz =");
    mpq_QSdump_array(p_mpq->lp->pIpiz, "pIpiz");
    QSlog("QSdelta_full_basis_status: pIdz =");
    mpq_QSdump_array(p_mpq->lp->pIdz, "pIdz");
  }

  mpq_ILLfct_check_pIdfeasible (p_mpq->lp, &fi, zero);
  mpq_ILLfct_set_status_values (p_mpq->lp, fi.pstatus, fi.dstatus,
                                           PHASEII,    PHASEI);
  if (p_mpq->lp->probstat.primal_feasible)
  {
    if (p_mpq->lp->probstat.dual_feasible)
    {
      *status = QS_LP_OPTIMAL;
    }
    else
    {
      *status = QS_LP_PRIMAL_FEASIBLE;
    }
  }
  else if (p_mpq->lp->probstat.primal_unbounded)
  {
    *status = QS_LP_UNBOUNDED;
  }
  else if (p_mpq->lp->probstat.primal_infeasible)
  {
    if (*status == QS_LP_INFEASIBLE)
      *simplexalgo = PRIMAL_SIMPLEX;  // More efficient than dual, if infeas
    *status = QS_LP_INFEASIBLE;
  }
  else if (p_mpq->lp->probstat.dual_feasible)
  {
      *status = QS_LP_DUAL_FEASIBLE;
  }
  else
    *status = QS_LP_UNSOLVED;
  EGtimerStop (&local_timer);
  if(!msg_lvl)
  {
    MESSAGE(0, "Performing Rational Basic Solve on %s, %s, check"
        " done in %lg seconds, PS %s %lg, DS %s %lg", p_mpq->name,
          *status == QS_LP_OPTIMAL         ? "RAT_optimal"
        : *status == QS_LP_UNBOUNDED       ? "RAT_unbounded"
        : *status == QS_LP_INFEASIBLE      ? "RAT_infeasible"
        : *status == QS_LP_PRIMAL_FEASIBLE ? "RAT_primal_feasible"
        : *status == QS_LP_DUAL_FEASIBLE   ? "RAT_dual_feasible"
                                           : "RAT_unsolved",
        local_timer.time,
          p_mpq->lp->basisstat.primal_feasible   ? "F"
        : p_mpq->lp->basisstat.primal_infeasible ? "I"
                                                 : "U",
          p_mpq->lp->basisstat.primal_feasible   ? mpq_get_d(p_mpq->lp->objval)
        : p_mpq->lp->basisstat.primal_infeasible ? mpq_get_d(p_mpq->lp->pinfeas)
                                                 : mpq_get_d(p_mpq->lp->objbound),
          p_mpq->lp->basisstat.dual_feasible   ? "F"
        : p_mpq->lp->basisstat.dual_infeasible ? "I"
                                               : "U",
          p_mpq->lp->basisstat.dual_feasible   ? mpq_get_d(p_mpq->lp->dobjval)
        : p_mpq->lp->basisstat.dual_infeasible ? mpq_get_d(p_mpq->lp->dinfeas)
                                               : mpq_get_d(p_mpq->lp->objbound));
  }
CLEANUP:
  mpq_EGlpNumClearVar (fi.totinfeas);
  return rval;
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
                        int *const simplexalgo,
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
  EGcallD(QSdelta_full_basis_status (p_mpq, status, basis, msg_lvl, simplexalgo));
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
    // FIXME: the new bound might not be *better*
    if (QS_LP_OPTIMAL == *status || QS_LP_PRIMAL_FEASIBLE == *status)
    {
      *have_primal = 1;
      EGcallD (copy_x (x, p_mpq));
      my_inner_prod (p_mpq->qslp->objsense == QS_MIN ? obj_up : obj_lo,
                            obj_coefs, x, p_mpq->qslp->nstruct);
    }
    if (QS_LP_OPTIMAL == *status || QS_LP_DUAL_FEASIBLE == *status)
    {
      *have_dual = 1;
      EGcallD (copy_y (y, p_mpq));
      my_inner_prod (p_mpq->qslp->objsense == QS_MIN ? obj_lo : obj_up,
                            rhs_coefs, y, p_mpq->qslp->nrows);
    }
    if (QS_LP_OPTIMAL == *status)
    {
      // Shortcut: can avoid difference check
      optimal_output (p_mpq);
      *judgement = 1;
      goto CLEANUP;
    }
    if (*have_primal && *have_dual)
    {
      // Could be optimal or delta-optimal
      mpq_sub (diff, obj_up, obj_lo);
      if (mpq_sgn (diff) == 0)
      {
        // Optimal
        optimal_output (p_mpq);
        *judgement = 1;
        goto CLEANUP;
      }
      else if (mpq_cmp (diff, delta) <= 0)
      {
        // Delta-optimal
        delta_optimal_output (p_mpq);
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
  mpq_EGlpNumInitVar (zero);  // Inits to zero

  if (!p_mpq || !p_mpq->qslp)
  {
    QSlog("LP data not available");
    rval = 1;
    goto CLEANUP;
  }
  int const msg_lvl = __QS_SB_VERB <= DEBUG ? 0: (1 - p_mpq->simplex_display) * 10000;
  *status = QS_LP_UNSOLVED;
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
    EGcallD (judge_basis (&judgement, 1, p_mpq, status, basis, msg_lvl, &simplexalgo,
                          delta, x, y, obj_lo, obj_up, obj_coefs, rhs_coefs,
                          &have_primal, &have_dual));
    if (judgement == 2)
    {
      MESSAGE(p_mpq->simplex_display ? 0 : __QS_SB_VERB, "double approximation"
              " failed, code %d, continuing in extended precision", rval);
      goto MPF_PRECISION;
    }
    else if (judgement)
      goto CLEANUP;
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
      EGcallD (judge_basis (&judgement, 0, p_mpq, status, basis, msg_lvl, &simplexalgo,
                            delta, x, y, obj_lo, obj_up, obj_coefs, rhs_coefs,
                            &have_primal, &have_dual));
      if (judgement)
        goto CLEANUP;
    }
    else
    {
      MESSAGE (msg_lvl, "Floating-point solver reports failure; trying next precision");
    }
  NEXT_PRECISION:
    mpf_QSfree_prob (p_mpf);
    p_mpf = 0;
  }

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
  mpq_QSfree_prob (p_mpq);

  EG_RETURN (rval);
}


/** @} */
