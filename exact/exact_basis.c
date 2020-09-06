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
/** @file
 * @ingroup Esolver */
/** @addtogroup Esolver */
/** @{ */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "exact_basis.h"
#include "dump.h"

#include <stdlib.h>
#include <string.h>

#include "logging-private.h"

#include "util.h"
#include "eg_timer.h"
#include "eg_exutil.h"
#include "except.h"

#include "basis_mpq.h"
#include "editor_dbl.h"
#include "editor_mpf.h"
#include "fct_mpq.h"
#include "lpdata_mpq.h"
#include "simplex_mpq.h"

/* ========================================================================= */
/** @brief get the status for a given basis in rational arithmetic, it should
 * also leave everything set to get primal/dual solutions when needed.
 * */
int QSexact_basis_status (mpq_QSdata * p_mpq,
                          int *status,
                          QSbasis * const basis,
                          const int msg_lvl,
                          int *const simplexalgo)
{
  int rval = 0,
  singular;
  mpq_feas_info fi;
  EGtimer_t local_timer;
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
  memset (&(p_mpq->lp->basisstat), 0, sizeof (mpq_lp_status_info));
  mpq_ILLfct_compute_piz (p_mpq->lp);
  mpq_ILLfct_compute_dz (p_mpq->lp);
  mpq_ILLfct_compute_xbz (p_mpq->lp);
  mpq_ILLfct_check_pfeasible (p_mpq->lp, &fi, mpq_zeroLpNum);
  mpq_ILLfct_check_dfeasible (p_mpq->lp, &fi, mpq_zeroLpNum);
  mpq_ILLfct_set_status_values (p_mpq->lp, fi.pstatus, fi.dstatus, PHASEII,
      PHASEII);
  if (p_mpq->lp->basisstat.optimal)
  {
    *status = QS_LP_OPTIMAL;
    EGcallD(mpq_QSgrab_cache (p_mpq, QS_LP_OPTIMAL));
  }
  else if (p_mpq->lp->basisstat.primal_infeasible
      || p_mpq->lp->basisstat.dual_unbounded)
  {
    if (*status == QS_LP_INFEASIBLE)
      *simplexalgo = PRIMAL_SIMPLEX;
    *status = QS_LP_INFEASIBLE;
    p_mpq->lp->final_phase = PRIMAL_PHASEI;
    p_mpq->lp->pIpiz = mpq_EGlpNumAllocArray (p_mpq->lp->nrows);
    mpq_ILLfct_compute_phaseI_piz (p_mpq->lp);
  }
  else if (p_mpq->lp->basisstat.primal_unbounded)
    *status = QS_LP_UNBOUNDED;
  else
    *status = QS_LP_UNSOLVED;
  EGtimerStop (&local_timer);
  if(!msg_lvl)
  {
    MESSAGE(0, "Performing Rational Basic Solve on %s, %s, check"
        " done in %lg seconds, PS %s %lg, DS %s %lg", p_mpq->name,
        (*status == QS_LP_OPTIMAL) ? "RAT_optimal" :
        ((*status == QS_LP_INFEASIBLE) ?  "RAT_infeasible" :
         ((*status == QS_LP_UNBOUNDED) ?  "RAT_unbounded" : "RAT_unsolved")),
        local_timer.time, p_mpq->lp->basisstat.primal_feasible ?
        "F":(p_mpq->lp->basisstat.primal_infeasible ? "I" : "U"),
        p_mpq->lp->basisstat.primal_feasible ?
        mpq_get_d(p_mpq->lp->objval) :
        (p_mpq->lp->basisstat.primal_infeasible ?
         mpq_get_d(p_mpq->lp->pinfeas) : mpq_get_d(p_mpq->lp->objbound)),
        p_mpq->lp->basisstat.dual_feasible ?
        "F":(p_mpq->lp->basisstat.dual_infeasible ? "I" : "U"),
        p_mpq->lp->basisstat.dual_feasible ? mpq_get_d(p_mpq->lp->dobjval)
        :(p_mpq->lp->basisstat.dual_infeasible ?
          mpq_get_d(p_mpq->lp->dinfeas) : mpq_get_d(p_mpq->lp->objbound)) );
  }
CLEANUP:
  mpq_EGlpNumClearVar (fi.totinfeas);
  return rval;
}

/* ========================================================================= */
/** @brief get the status for a given basis in rational arithmetic, it should
 * also leave everything set to get primal/dual solutions when needed.
 * */
int QSdelta_basis_status (mpq_QSdata * p_mpq,
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
    QSlog("QSdelta_basis_status: xnbz =");
    mpq_QSdump_xnbz(p_mpq);
    QSlog("QSdelta_basis_status: xbz =");
    mpq_QSdump_xbz(p_mpq);
    QSlog("QSdelta_basis_status: bfeas =");
    mpq_QSdump_bfeas(p_mpq);
    QSlog("QSdelta_basis_status: pIpiz =");
    mpq_QSdump_array(p_mpq->lp->pIpiz, "pIpiz");
    QSlog("QSdelta_basis_status: pIdz =");
    mpq_QSdump_array(p_mpq->lp->pIdz, "pIdz");
  }

  mpq_ILLfct_check_pIdfeasible (p_mpq->lp, &fi, zero);
  mpq_ILLfct_set_status_values (p_mpq->lp, fi.pstatus, fi.dstatus,
                                           PHASEII,    PHASEI);
  if (p_mpq->lp->probstat.primal_feasible
   || p_mpq->lp->probstat.primal_unbounded)
    *status = QS_LP_FEASIBLE;
  else if (p_mpq->lp->probstat.primal_infeasible)
  {
    if (*status == QS_LP_INFEASIBLE)
      *simplexalgo = PRIMAL_SIMPLEX;  // More efficient than dual, if infeas
    *status = QS_LP_INFEASIBLE;
  }
  else
    *status = QS_LP_UNSOLVED;
  EGtimerStop (&local_timer);
  if(!msg_lvl)
  {
    MESSAGE(0, "Performing Rational Basic Solve on %s, %s, check"
        " done in %lg seconds, PS %s %lg, DS %s %lg", p_mpq->name,
          *status == QS_LP_FEASIBLE   ? "RAT_feasible"
        : *status == QS_LP_INFEASIBLE ? "RAT_infeasible"
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

/* ========================================================================= */
/** @} */
/* end of exact_basis.c */

