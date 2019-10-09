/** @file
 *  @ingroup DeltaSolver */
/** @addtogroup DeltaSolver
 *  @{ */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "delta.h"

#include "exact.h"
#include "except.h"
#include "editor_dbl.h"
#include "editor_mpf.h"
#include "eg_macros.h"

/* ========================================================================= */
/** @brief Used as separator while printing output to the screen (controlled by
 * enabling simplex_display in the mpq_QSdata */
/* ========================================================================= */
static const char __sp[81] =
  "================================================================================";

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
    QSlog("Problem Is Infeasible");
  }
  if (y)
  {
    unsigned sz = __EGlpNumArraySize (y_mpq);
    while (sz--)
      mpq_set (y[sz], y_mpq[sz]);
  }
}

/* ========================================================================= */
/** @brief print into screen (if enable) a message indicating that we have
 * successfully solved the problem at optimality, and save (if x and y are non
 * NULL respectivelly) the optimal primal/dual solution provided in x_mpq and
 * y_mpq. 
 * @param p_mpq the problem data.
 * @param x where to store the optimal primal solution (if not null).
 * @param y where to store the optimal dual solution (if not null).
 * @param x_mpq  the optimal primal solution.
 * @param y_mpq  the optimal dual solution.
 * */
/* ========================================================================= */
static void optimal_output (mpq_QSdata * p_mpq,
                            mpq_t * const x,
                            mpq_t * const y,
                            mpq_t * x_mpq,
                            mpq_t * y_mpq)
{
  if (p_mpq->simplex_display)
  {
    QSlog("Problem Solved Exactly");
  }
  if (y)
  {
    unsigned sz = __EGlpNumArraySize (y_mpq);
    while (sz--)
      mpq_set (y[sz], y_mpq[sz]);
  }
  if (x)
  {
    unsigned sz = __EGlpNumArraySize (x_mpq);
    while (sz--)
      mpq_set (x[sz], x_mpq[sz]);
  }
}

/* ========================================================================= */
/** @brief Given an mpq_QSdata problem, solve the corresponding
 * delta-feasibility problem exactly.
 * @param p_mpq problem for which to determine delta-feasibility exactly.
 * @param delta the delta to use for determining delta-feasibility; the maximum
 * perturbation of RHS/bounds required to make a delta-feasible solution
 * feasible.
 * @param x if not null, we store here a delta-feasible solution to the problem
 * (if delta-feasibility established).
 * @param y if not null, we store here a certificate of infeasibility for the
 * problem (if infeasibility established).
 * @param ebasis if not null, use the given basis to start the iteration of
 * simplex, and store here the final basis (where applicable).
 * @param simplexalgo whether to use primal or dual simplex while solving the
 * delta-feasibility problem.
 * @param status pointer to the integer where we will return the status of the
 * problem, either delta-feasible or infeasible (we could also return time
 * out).
 * @return zero on success, non-zero otherwise. */
int QSdelta_solver (mpq_QSdata * p_mpq,
                    mpq_t const delta,
                    mpq_t * const x,
                    mpq_t * const y,
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
  mpq_t *x_mpq = 0,
   *y_mpq = 0;
  mpf_t *x_mpf = 0,
   *y_mpf = 0;
  int const msg_lvl = __QS_SB_VERB <= DEBUG ? 0: (1 - p_mpq->simplex_display) * 10000;
  *status = 0;
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
  if(__QS_SB_VERB <= DEBUG) p_dbl->simplex_display = 1;
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
  switch (*status)
  {
  case QS_LP_OPTIMAL:
    x_dbl = dbl_EGlpNumAllocArray (p_dbl->qslp->ncols);
    y_dbl = dbl_EGlpNumAllocArray (p_dbl->qslp->nrows);
    EGcallD(dbl_QSget_x_array (p_dbl, x_dbl));
    EGcallD(dbl_QSget_pi_array (p_dbl, y_dbl));
    x_mpq = QScopy_array_dbl_mpq (x_dbl);
    y_mpq = QScopy_array_dbl_mpq (y_dbl);
    dbl_EGlpNumFreeArray (x_dbl);
    dbl_EGlpNumFreeArray (y_dbl);
    basis = dbl_QSget_basis (p_dbl);
    if (QSexact_optimal_test (p_mpq, x_mpq, y_mpq, basis))
    {
      optimal_output (p_mpq, x, y, x_mpq, y_mpq);
      goto CLEANUP;
    }
    else
    {
      EGcallD(QSexact_basis_status (p_mpq, status, basis, msg_lvl, &simplexalgo));
      if (*status == QS_LP_OPTIMAL)
      {
        if(!msg_lvl)
        {
          MESSAGE(0,"Retesting solution");
        }
        EGcallD(mpq_QSget_x_array (p_mpq, x_mpq));
        EGcallD(mpq_QSget_pi_array (p_mpq, y_mpq));
        if (QSexact_optimal_test (p_mpq, x_mpq, y_mpq, basis))
        {
          optimal_output (p_mpq, x, y, x_mpq, y_mpq);
          goto CLEANUP;
        }
        else
        {
          last_status = *status = QS_LP_UNSOLVED;
        }
      }
      else
      {
        if(!msg_lvl)
        {
          MESSAGE(0,"Status is not optimal, but %d", *status);
        }
      }
    }
    mpq_EGlpNumFreeArray (x_mpq);
    mpq_EGlpNumFreeArray (y_mpq);
    break;
  case QS_LP_INFEASIBLE:
    y_dbl = dbl_EGlpNumAllocArray (p_dbl->qslp->nrows);
    if (dbl_QSget_infeas_array (p_dbl, y_dbl))
    {
      MESSAGE(p_mpq->simplex_display ? 0 : __QS_SB_VERB, "double approximation"
              " failed, code %d, continuing in extended precision\n", rval);
      goto MPF_PRECISION;
    }
    y_mpq = QScopy_array_dbl_mpq (y_dbl);
    dbl_EGlpNumFreeArray (y_dbl);
    if (QSexact_infeasible_test (p_mpq, y_mpq))
    {
      infeasible_output (p_mpq, y, y_mpq);
      goto CLEANUP;
    }
    else
    {
      MESSAGE (msg_lvl, "Retesting solution in exact arithmetic");
      basis = dbl_QSget_basis (p_dbl);
      EGcallD(QSexact_basis_status (p_mpq, status, basis, msg_lvl, &simplexalgo));
      if (*status == QS_LP_INFEASIBLE)
      {
        mpq_EGlpNumFreeArray (y_mpq);
        y_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->nrows);
        EGcallD(mpq_QSget_infeas_array (p_mpq, y_mpq));
        if (QSexact_infeasible_test (p_mpq, y_mpq))
        {
          infeasible_output (p_mpq, y, y_mpq);
          goto CLEANUP;
        }
        else
        {
          last_status = *status = QS_LP_UNSOLVED;
        }
      }
    }
    mpq_EGlpNumFreeArray (y_mpq);
    break;
  case QS_LP_UNBOUNDED:
    MESSAGE(p_mpq->simplex_display ? 0 : __QS_SB_VERB, "%s\n\tUnbounded "
            "Problem found, not implemented to deal with this\n%s\n",__sp,__sp);
    break;
  case QS_LP_OBJ_LIMIT:
    rval=1;
    IFMESSAGE(p_mpq->simplex_display,"Objective limit reached (in floating point) ending now");
    goto CLEANUP;
    break;
  default:
    IFMESSAGE(p_mpq->simplex_display,"Re-trying inextended precision");
    break;
  }
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
    if(__QS_SB_VERB <= DEBUG) p_mpf->simplex_display = 1;
    simplexalgo = PRIMAL_SIMPLEX;
    if(!last_iter) last_status = QS_LP_UNSOLVED;
    if(last_status == QS_LP_OPTIMAL || last_status == QS_LP_INFEASIBLE)
    {
      if (p_mpq->simplex_display || DEBUG >= __QS_SB_VERB)
      {
        QSlog("Re-using previous basis");
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
        QSlog("Not-using previous basis");
      }
    }
    if (mpf_ILLeditor_solve (p_mpf, simplexalgo))
    {
      if (p_mpq->simplex_display || DEBUG >= __QS_SB_VERB)
      {
        QSlog("mpf_%u precision falied, error code %d, continuing with "
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
    switch (*status)
    {
    case QS_LP_OPTIMAL:
      basis = mpf_QSget_basis (p_mpf);
      x_mpf = mpf_EGlpNumAllocArray (p_mpf->qslp->ncols);
      y_mpf = mpf_EGlpNumAllocArray (p_mpf->qslp->nrows);
      EGcallD(mpf_QSget_x_array (p_mpf, x_mpf));
      EGcallD(mpf_QSget_pi_array (p_mpf, y_mpf));
      x_mpq = QScopy_array_mpf_mpq (x_mpf);
      y_mpq = QScopy_array_mpf_mpq (y_mpf);
      mpf_EGlpNumFreeArray (x_mpf);
      mpf_EGlpNumFreeArray (y_mpf);
      if (QSexact_optimal_test (p_mpq, x_mpq, y_mpq, basis))
      {
        optimal_output (p_mpq, x, y, x_mpq, y_mpq);
        goto CLEANUP;
      }
      else
      {
        EGcallD(QSexact_basis_status (p_mpq, status, basis, msg_lvl, &simplexalgo));
        if (*status == QS_LP_OPTIMAL)
        {
          MESSAGE (msg_lvl, "Retesting solution");
          EGcallD(mpq_QSget_x_array (p_mpq, x_mpq));
          EGcallD(mpq_QSget_pi_array (p_mpq, y_mpq));
          if (QSexact_optimal_test (p_mpq, x_mpq, y_mpq, basis))
          {
            optimal_output (p_mpq, x, y, x_mpq, y_mpq);
            goto CLEANUP;
          }
          else
          {
            last_status = *status = QS_LP_UNSOLVED;
          }
        }
        else
          MESSAGE (msg_lvl, "Status is not optimal, but %d", *status);
      }
      mpq_EGlpNumFreeArray (x_mpq);
      mpq_EGlpNumFreeArray (y_mpq);
      break;
    case QS_LP_INFEASIBLE:
      y_mpf = mpf_EGlpNumAllocArray (p_mpf->qslp->nrows);
      EGcallD(mpf_QSget_infeas_array (p_mpf, y_mpf));
      y_mpq = QScopy_array_mpf_mpq (y_mpf);
      mpf_EGlpNumFreeArray (y_mpf);
      if (QSexact_infeasible_test (p_mpq, y_mpq))
      {
        infeasible_output (p_mpq, y, y_mpq);
        goto CLEANUP;
      }
      else
      {
        MESSAGE (msg_lvl, "Retesting solution in exact arithmetic");
        basis = mpf_QSget_basis (p_mpf);
        EGcallD(QSexact_basis_status (p_mpq, status, basis, msg_lvl, &simplexalgo));
        if (*status == QS_LP_INFEASIBLE)
        {
          mpq_EGlpNumFreeArray (y_mpq);
          y_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->nrows);
          EGcallD(mpq_QSget_infeas_array (p_mpq, y_mpq));
          if (QSexact_infeasible_test (p_mpq, y_mpq))
          {
            infeasible_output (p_mpq, y, y_mpq);
            goto CLEANUP;
          }
          else
          {
            last_status = *status = QS_LP_UNSOLVED;
          }
        }
      }
      mpq_EGlpNumFreeArray (y_mpq);
      break;
      break;
    case QS_LP_OBJ_LIMIT:
      rval=1;
      IFMESSAGE(p_mpq->simplex_display,"Objective limit reached (in floating point) ending now");
      goto CLEANUP;
      break;
    case QS_LP_UNBOUNDED:
    default:
      MESSAGE(__QS_SB_VERB,"Re-trying inextended precision");
      break;
    }
  NEXT_PRECISION:
    mpf_QSfree_prob (p_mpf);
    p_mpf = 0;
  }
  /* ending */
CLEANUP:
  dbl_EGlpNumFreeArray (x_dbl);
  dbl_EGlpNumFreeArray (y_dbl);
  mpq_EGlpNumFreeArray (x_mpq);
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
  return rval;
}

/** @} */
