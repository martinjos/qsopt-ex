/** @file
 *  @ingroup DeltaSolver */
/** @addtogroup DeltaSolver
 *  @{ */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "delta.h"
#include "exact_basis.h"

#include "qstruct_dbl.h"
#include "qstruct_mpf.h"

#include "basis_mpq.h"
#include "except.h"
#include "editor_dbl.h"
#include "editor_mpf.h"
#include "eg_macros.h"
#include "eg_timer.h"
#include "fct_mpq.h"
#include "qsopt_mpf.h"
#include "qsopt_dbl.h"
#include "simplex_mpq.h"
#include "dump.h"

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
                               mpq_t * y_mpq)
{
  if (p_mpq->simplex_display)
  {
    QSlog("Problem is infeasible");
  }
  if (y)
  {
    unsigned sz = __EGlpNumArraySize (y_mpq);
    while (sz--)
      mpq_set (y[sz], y_mpq[sz]);
  }
}

int QSdelta_copy_x (mpq_t * const x, const mpq_QSdata * p_mpq)
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
    EGcallD(QSdelta_copy_x (x, p_mpq));

CLEANUP:

  EG_RETURN (rval);
}

/* ========================================================================= */
/** Check for and report delta-feasibility of the basis
 * @param p_mpq the problem data.
 * @param delta the maximum infeasibility of a delta-feasible solution; updated
 * with the actual infeasibility if delta-feasibility determined.
 * @param status where to store the new status (feasible or delta-feasible),
 * if applicable
 * @param x where to store a delta-feasible primal solution (if not null).
 * */
/* ========================================================================= */
static int check_delta_feas (mpq_QSdata const * p_mpq,
                             mpq_t delta,
                             int *status,
                             mpq_t * const x,
                             delta_callback_t delta_callback,
                             mpq_t last_infeas,
                             void *callback_data)
{
  int i, col;
  mpq_t infeas, err1, err2;
  int rval = 0;
  mpq_lpinfo* lp = p_mpq->lp;
  mpq_ILLlpdata* qslp = p_mpq->qslp;

  *status = QS_LP_UNSOLVED;

  mpq_EGlpNumInitVar (infeas);
  mpq_EGlpNumInitVar (err1);
  mpq_EGlpNumInitVar (err2);
  mpq_EGlpNumZero (infeas);

  for (i = 0; i < lp->nrows; i++)
  {
    col = lp->baz[i];
    mpq_EGlpNumCopyDiff (err1, lp->xbz[i], lp->uz[col]);
    mpq_EGlpNumCopyDiff (err2, lp->lz[col], lp->xbz[i]);
    if (mpq_EGlpNumIsLess (mpq_zeroLpNum, err1)
        && mpq_EGlpNumIsNeq (lp->uz[col], mpq_INFTY, mpq_zeroLpNum))
    {
      if (mpq_EGlpNumIsLess (infeas, err1))
        mpq_EGlpNumCopy (infeas, err1);
      WARNINGL (QSE_WLVL, mpq_EGlpNumIsLess (mpq_INFTY, err1),
               "This is impossible: lu = %15lg xbz = %15lg" " mpq_INFTY = %15lg",
               mpq_EGlpNumToLf (lp->uz[col]), mpq_EGlpNumToLf (lp->xbz[i]),
               mpq_EGlpNumToLf (mpq_INFTY));
    }
    else if (mpq_EGlpNumIsLess (mpq_zeroLpNum, err2)
             && mpq_EGlpNumIsNeq (lp->lz[col], mpq_NINFTY, mpq_zeroLpNum))
    {
      if (mpq_EGlpNumIsLess (infeas, err2))
        mpq_EGlpNumAddTo (infeas, err2);
      WARNINGL (QSE_WLVL, mpq_EGlpNumIsLess (mpq_INFTY, err2),
               "This is impossible: lz = %15lg xbz = %15lg" " mpq_NINFTY = %15lg",
               mpq_EGlpNumToLf (lp->lz[col]), mpq_EGlpNumToLf (lp->xbz[i]),
               mpq_EGlpNumToLf (mpq_NINFTY));
    }
  }

  if (mpq_EGlpNumIsLessZero (infeas))
  {
    QSlog("Negative infeasibility (impossible): %lf %la",
                mpq_EGlpNumToLf (infeas), mpq_EGlpNumToLf (infeas));
  }

  if (!mpq_EGlpNumIsNeqqZero (infeas))
  {
    // feasible
    if (p_mpq->simplex_display)
    {
      QSlog("Problem is feasible");
    }
    mpq_EGlpNumCopy (delta, infeas);
    *status = QS_LP_FEASIBLE;
  }
  else if (!mpq_EGlpNumIsLess (delta, infeas))
  {
    // delta-feasible
    if (p_mpq->simplex_display)
    {
      QSlog("Problem is delta-feasible with delta = %lf",
            mpq_EGlpNumToLf (infeas));
    }
    mpq_EGlpNumCopy (delta, infeas);
    *status = QS_LP_DELTA_FEASIBLE;
  }
  else if (NULL != delta_callback)
  {
    if (mpq_sgn (last_infeas) == 0 || mpq_cmp (infeas, last_infeas) < 0)
    {
      mpq_set (last_infeas, infeas);
      delta_callback(p_mpq, x, infeas, delta, callback_data);
    }
  }

  if (x && (QS_LP_FEASIBLE == *status || QS_LP_DELTA_FEASIBLE == *status))
    EGcallD(QSdelta_copy_x (x, p_mpq));

  if (QS_LP_FEASIBLE != *status && QS_LP_DELTA_FEASIBLE != *status)
  {
    IFMESSAGE(p_mpq->simplex_display, "Failed to make final conclusion on basis of solver result");
  }

CLEANUP:

  mpq_EGlpNumClearVar (infeas);
  mpq_EGlpNumClearVar (err1);
  mpq_EGlpNumClearVar (err2);

  EG_RETURN (rval);
}

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
 * @param delta_callback if not null, will be called if a delta-satisfying
 * result is found for some value greater than delta.
 * @param callback_data additional parameter to be passed to delta_callback.
 * @return zero on success, non-zero otherwise. */
int QSdelta_solver (mpq_QSdata * p_orig,
                    mpq_t delta,
                    mpq_t * const x,
                    mpq_t * const y,
                    QSbasis * const ebasis,
                    int simplexalgo,
                    int *status,
                    delta_callback_t delta_callback,
                    void *callback_data)
{
  /* local variables */
  int last_status = 0, last_iter = 0;
  QSbasis *basis = 0;
  unsigned precision = EGLPNUM_PRECISION;
  int rval = 0,
    it = QS_EXACT_MAX_ITER;
  mpq_QSdata *p_mpq = 0;
  dbl_QSdata *p_dbl = 0;
  mpf_QSdata *p_mpf = 0;
  double *x_dbl = 0,
   *y_dbl = 0;
  mpq_t *y_mpq = 0;
  mpf_t *x_mpf = 0,
   *y_mpf = 0;
  mpq_t last_infeas;
  mpq_init (last_infeas);
  int const msg_lvl = __QS_SB_VERB <= DEBUG ? 0: (1 - p_orig->simplex_display) * 10000;
  *status = QS_LP_UNSOLVED;
  p_mpq = mpq_QScopy_prob (p_orig, "mpq_feas_problem");
  /* set the objective function to zero (in the copy) */
  EGcallD (mpq_QSclear_obj (p_mpq));
  /* save the problem if we are really debugging */
  if(DEBUG >= __QS_SB_VERB)
  {
    EGcallD(mpq_QSwrite_prob(p_mpq, "qsxprob.lp","LP"));
  }
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
    MESSAGE (msg_lvl, "Basis hash is 0x%016lX", QSexact_basis_hash(basis));
    EGcallD(QSdelta_basis_status (p_mpq, status, basis, msg_lvl, &simplexalgo));
    if (QS_LP_INFEASIBLE == *status)
    {
      y_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->nrows);
      if (mpq_QSget_infeas_array (p_mpq, y_mpq))
      {
        MESSAGE(p_mpq->simplex_display ? 0 : __QS_SB_VERB, "double approximation"
                " failed, code %d, continuing in extended precision", rval);
        goto MPF_PRECISION;
      }
      infeasible_output (p_mpq, y, y_mpq);
      goto CLEANUP;
    }
    else if (QS_LP_FEASIBLE == *status)
    {
      EGcallD(feasible_output (p_mpq, x));
      mpq_EGlpNumCopy (delta, mpq_zeroLpNum);
      goto CLEANUP;
    }
    /* check for delta-feasibility */
    EGcallD(check_delta_feas (p_mpq, delta, status, x, delta_callback, last_infeas, callback_data));
    if (QS_LP_FEASIBLE == *status || QS_LP_DELTA_FEASIBLE == *status)
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
      MESSAGE (msg_lvl, "Basis hash is 0x%016lX", QSexact_basis_hash(basis));
      EGcallD(QSdelta_basis_status (p_mpq, status, basis, msg_lvl, &simplexalgo));
      if (QS_LP_INFEASIBLE == *status)
      {
        mpq_EGlpNumFreeArray (y_mpq);
        y_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->nrows);
        EGcallD(mpq_QSget_infeas_array (p_mpq, y_mpq));
        infeasible_output (p_mpq, y, y_mpq);
        goto CLEANUP;
      }
      else if (QS_LP_FEASIBLE == *status)
      {
        EGcallD(feasible_output (p_mpq, x));
        mpq_EGlpNumCopy (delta, mpq_zeroLpNum);
        goto CLEANUP;
      }
      /* check for delta-feasibility */
      EGcallD(check_delta_feas (p_mpq, delta, status, x, delta_callback, last_infeas, callback_data));
      if (QS_LP_FEASIBLE == *status || QS_LP_DELTA_FEASIBLE == *status)
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
  /* ending */
CLEANUP:
  mpq_clear (last_infeas);
  dbl_EGlpNumFreeArray (x_dbl);
  dbl_EGlpNumFreeArray (y_dbl);
  mpq_EGlpNumFreeArray (y_mpq);
  mpf_EGlpNumFreeArray (x_mpf);
  mpf_EGlpNumFreeArray (y_mpf);
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
