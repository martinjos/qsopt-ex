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
/** @file
 * @ingroup Esolver */
/** @addtogroup Esolver */
/** @{ */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "exact_delta.h"
#include "exact_delta_g_mpq.h"
#include "exact_delta_g_mpf.h"
#include "exact_delta_g_dbl.h"

#include "logging-private.h"

#include "except.h"

#include "editor_dbl.h"
#include "editor_mpf.h"
#include "lpdata_mpq.h"
#include "fct_mpq.h"

/* ========================================================================= */
int QSexact_delta_optimal_test (mpq_QSdata * p,
																mpq_t * p_sol,
																mpq_t * d_sol,
																QSbasis * basis,
																mpq_t const delta)
{
	/* local variables */
	register int i,
	  j;
	mpq_ILLlpdata *qslp = p->lp->O;
	int *rowmap = qslp->rowmap,
	 *structmap = qslp->structmap,
	  col;
	mpq_t *rhs_copy = 0;
	mpq_t *dz = 0;
	int objsense = (qslp->objsense == QS_MIN) ? 1 : -1;
	int const msg_lvl = __QS_SB_VERB <= DEBUG ? 0 : 100000 * (1 - p->simplex_display);
	int rval = QS_EXACT_UNKNOWN;  /* safe fallback value */
	mpq_t num1,
	  num2,
	  num3,
	  d_obj;
	mpq_init (num1);
	mpq_init (num2);
	mpq_init (num3);
	mpq_init (d_obj);
	mpq_set_ui (d_obj, 0UL, 1UL);

	if(!msg_lvl)
	{
		unsigned sz;
		QSlog("QSexact_delta_optimal_test: p_sol =");
		sz = __EGlpNumArraySize (p_sol);
		for (int i = 0; i < sz; ++i)
		{
			QSlog("%d: %g", i, mpq_get_d (p_sol[i]));
		}
		QSlog("QSexact_delta_optimal_test: d_sol =");
		sz = __EGlpNumArraySize (d_sol);
		for (int i = 0; i < sz; ++i)
		{
			QSlog("%d: %g", i, mpq_get_d (d_sol[i]));
		}
	}

	/* now check if the given basis is the optimal basis */
	if (mpq_QSload_basis (p, basis))
	{
		rval = QS_EXACT_UNKNOWN;
		MESSAGE (msg_lvl, "QSload_basis failed");
		goto CLEANUP;
	}
	for (i = basis->nstruct; i--;)
	{
		/* check that the upper and lower bound define a non-empty space */
		if (mpq_cmp (qslp->lower[structmap[i]], qslp->upper[structmap[i]]) > 0)
		{
			rval = QS_EXACT_UNSAT;
			if(!msg_lvl)
			{
				MESSAGE(0, "variable %s has empty feasible range [%lg,%lg]",
								 qslp->colnames[i], mpq_EGlpNumToLf(qslp->lower[structmap[i]]),
								 mpq_EGlpNumToLf(qslp->upper[structmap[i]]));
			}
			goto CLEANUP;
		}
		/* set the variable to its apropiate values, depending its status */
		if (0 != mpq_sgn (qslp->obj[structmap[i]]))
		{
			/* artificial variable, added by my solver - set value to zero */
			mpq_set_ui(p_sol[i], 0UL, 1UL);
			if (0 != mpq_sgn(qslp->lower[structmap[i]]))
			{
				rval = QS_EXACT_UNKNOWN;
				MESSAGE(0, "ERROR IN INPUT: variable %s has non-zero objective"
				           " coefficient, and its lower bound is not zero",
								 qslp->colnames[i]);
				goto CLEANUP;
			}
		}
		else
		{
			switch (basis->cstat[i])
			{
			case QS_COL_BSTAT_FREE:
			case QS_COL_BSTAT_BASIC:
				if (mpq_cmp (p_sol[i], qslp->upper[structmap[i]]) > 0)
					mpq_set (p_sol[i], qslp->upper[structmap[i]]);
				else if (mpq_cmp (p_sol[i], qslp->lower[structmap[i]]) < 0)
					mpq_set (p_sol[i], qslp->lower[structmap[i]]);
				break;
			case QS_COL_BSTAT_UPPER:
				mpq_set (p_sol[i], qslp->upper[structmap[i]]);
				break;
			case QS_COL_BSTAT_LOWER:
				mpq_set (p_sol[i], qslp->lower[structmap[i]]);
				break;
			default:
				rval = QS_EXACT_UNKNOWN;
				MESSAGE (msg_lvl, "Unknown Variable basic status %d, for variable "
								 "(%s,%d)", basis->cstat[i], qslp->colnames[i], i);
				goto CLEANUP;
				break;
			}
		}
	}
	for (i = basis->nrows; i--;)
	{
		/* check that the upper and lower bound define a non-empty space */
		if (mpq_cmp (qslp->lower[rowmap[i]], qslp->upper[rowmap[i]]) > 0)
		{
			rval = QS_EXACT_UNSAT;
			if(!msg_lvl)
			{
				MESSAGE(0, "constraint %s logical has empty feasible range "
								 "[%lg,%lg]", qslp->rownames[i], 
								 mpq_EGlpNumToLf(qslp->lower[rowmap[i]]),
								 mpq_EGlpNumToLf(qslp->upper[rowmap[i]]));
			}
			goto CLEANUP;
		}
		/* no need to set the value here, as it would get overwritten later */
	}

	/* compute the actual RHS */
	rhs_copy = mpq_EGlpNumAllocArray (qslp->nrows);
	for (i = qslp->nstruct; i--;)
	{
		if (!mpq_equal (p_sol[i], mpq_zeroLpNum))
		{
			mpq_t* arr = qslp->A.matval + qslp->A.matbeg[structmap[i]];
			int* iarr = qslp->A.matind + qslp->A.matbeg[structmap[i]];
			for (j = qslp->A.matcnt[structmap[i]]; j--;)
			{
				mpq_mul (num1, arr[j], p_sol[i]);
				mpq_add (rhs_copy[iarr[j]], rhs_copy[iarr[j]], num1);
			}
		}
	}

	rval = QS_EXACT_SAT;  /* provisionally SAT */

	/* now replace the row slack, and check for delta-sat/sat */
	for (i = qslp->nrows; i--;)
	{
		mpq_mul (num1, qslp->rhs[i], d_sol[i]);
		mpq_add (d_obj, d_obj, num1);
		mpq_sub (num2, qslp->rhs[i], rhs_copy[i]);  /* num2 = b_i - A_i.x */
		/* clamp to "infinities" */
		if (mpq_cmp (num2, mpq_INFTY) > 0)
		{
			mpq_set(num2, mpq_INFTY);
		}
		else if (mpq_cmp (num2, mpq_NINFTY) < 0)
		{
			mpq_set(num2, mpq_NINFTY);
		}
		EXIT (qslp->A.matcnt[rowmap[i]] != 1, "Imposible!");
		/* divide by its own coefficient to get value of row slack */
		mpq_div (num2, num2,
						 qslp->A.matval[qslp->A.matbeg[rowmap[i]]]);
		/* always replace row slacks, even if non-basic */
		mpq_set (p_sol[qslp->nstruct + i], num2);
		/* now we check the bounds on the logical variables */
		if (QS_EXACT_UNKNOWN != rval)
		{
			if (mpq_cmp (num2, qslp->lower[rowmap[i]]) < 0)
			{
				mpq_add(num2, num2, delta);
				if (mpq_cmp (num2, qslp->lower[rowmap[i]]) < 0)
					rval = QS_EXACT_UNKNOWN;  /* may be unsat */
				else
					rval = QS_EXACT_DELTA_SAT;  /* provisionally delta-sat */
			}
			else if (mpq_cmp (num2, qslp->upper[rowmap[i]]) > 0)
			{
				mpq_sub(num2, num2, delta);
				if (mpq_cmp (num2, qslp->upper[rowmap[i]]) > 0)
					rval = QS_EXACT_UNKNOWN;  /* may be unsat */
				else
					rval = QS_EXACT_DELTA_SAT;  /* provisionally delta-sat */
			}
		}
	}

	if (QS_EXACT_SAT == rval || QS_EXACT_DELTA_SAT == rval)
	{
		/* we already have proof */
		if(!msg_lvl)
		{
			MESSAGE(0, "Problem is proven %ssatisfiable",
							QS_EXACT_DELTA_SAT==rval ? "delta-" : "");
		}
		goto CLEANUP;
	}
	else
	{
		if(!msg_lvl)
		{
			MESSAGE(0, "Failed to prove problem delta-satisfiable");
		}
	}

	/* We have now determined that nothing can be learnt from the primal
	 * solution. If the dual vector is feasible and the dual objective value is
	 * positive, we can conclude unsat; otherwise, the status is unknown. */

	/* finish computing the dual objective function */
	/* compute the upper and lower bound dual variables, note that dl is the dual
	 * of the lower bounds, and du the dual of the upper bound, dl >= 0 and du <=
	 * 0 and A^t y + Idl + Idu = c, and the dual objective value is 
	 * max y*b + l*dl + u*du, we colapse both vector dl and du into dz, note that
	 * if we are maximizing, then dl <= 0 and du >=0 */
	dz = mpq_EGlpNumAllocArray (qslp->ncols);
	for (i = qslp->nstruct; i--;)
	{
		col = structmap[i];
		mpq_t* arr = qslp->A.matval + qslp->A.matbeg[col];
		int* iarr = qslp->A.matind + qslp->A.matbeg[col];
		mpq_set (num1, qslp->obj[col]);
		for (j = qslp->A.matcnt[col]; j--;)
		{
			mpq_mul (num2, arr[j], d_sol[iarr[j]]);
			mpq_sub (num1, num1, num2);
		}
		mpq_set (dz[col], num1);
		/* objective update */
		if (objsense * mpq_cmp_ui (dz[col], 0UL, 1UL) > 0)
		{
			mpq_mul (num3, dz[col], qslp->lower[col]);
			mpq_add (d_obj, d_obj, num3);
		}
		else
		{
			mpq_mul (num3, dz[col], qslp->upper[col]);
			mpq_add (d_obj, d_obj, num3);
		}
		/* dual feasibility */
		int dual_sense = objsense * mpq_sgn (dz[col]);
		if (dual_sense > 0 && mpq_equal(qslp->lower[col], mpq_NINFTY))
		{
			/* dual infeasible */
			rval = QS_EXACT_UNKNOWN;
			if(!msg_lvl)
			{
				MESSAGE(0, "solution is dual infeasible for variable %s"
									 " (lower bound infinite)", qslp->colnames[i]);
			}
			goto CLEANUP;
		}
		else if (dual_sense < 0 && mpq_equal(qslp->upper[col], mpq_INFTY))
		{
			/* dual infeasible */
			rval = QS_EXACT_UNKNOWN;
			if(!msg_lvl)
			{
				MESSAGE(0, "solution is dual infeasible for variable %s"
									 " (upper bound infinite)", qslp->colnames[i]);
			}
			goto CLEANUP;
		}
	}
	/* finished problem variables, now do the same for the logical
	 * variables (i.e. row slacks) */
	for (i = qslp->nrows; i--;)
	{
		col = rowmap[i];
		WARNING (mpq_cmp (qslp->obj[col], mpq_zeroLpNum),
						 "logical variable %s with non-zero objective function %lf",
						 qslp->rownames[i], mpq_get_d (qslp->obj[col]));
		mpq_t* arr = qslp->A.matval + qslp->A.matbeg[col];
		int* iarr = qslp->A.matind + qslp->A.matbeg[col];
		mpq_set (num1, qslp->obj[col]);
		for (j = qslp->A.matcnt[col]; j--;)
		{
			mpq_mul (num2, arr[j], d_sol[iarr[j]]);
			mpq_sub (num1, num1, num2);
		}
		mpq_set (dz[col], num1);
		/* objective update */
		if (objsense * mpq_cmp_ui (dz[col], 0UL, 1UL) > 0)
		{
			mpq_mul (num3, dz[col], qslp->lower[col]);
			mpq_add (d_obj, d_obj, num3);
		}
		else
		{
			mpq_mul (num3, dz[col], qslp->upper[col]);
			mpq_add (d_obj, d_obj, num3);
		}
		/* dual feasibility */
		int dual_sense = objsense * mpq_sgn (dz[col]);
		if (dual_sense > 0 && mpq_equal(qslp->lower[col], mpq_NINFTY))
		{
			/* dual infeasible */
			rval = QS_EXACT_UNKNOWN;
			if(!msg_lvl)
			{
				MESSAGE(0, "solution is dual infeasible for logical variable %s"
									 " (lower bound infinite)", qslp->rownames[i]);
			}
			goto CLEANUP;
		}
		else if (dual_sense < 0 && mpq_equal(qslp->upper[col], mpq_INFTY))
		{
			/* dual infeasible */
			rval = QS_EXACT_UNKNOWN;
			if(!msg_lvl)
			{
				MESSAGE(0, "solution is dual infeasible for logical variable %s"
									 " (upper bound infinite)", qslp->rownames[i]);
			}
			goto CLEANUP;
		}
	}

	/* d_obj now provides a rigorous lower bound on the optimal objective value.
	 * If this is positive, the original problem is unsatisfiable. */
	if (mpq_cmp_ui (d_obj, 0UL, 1UL) > 0)
	{
		if(!msg_lvl)
		{
			MESSAGE(0, "Problem is proven unsatisfiable (dual objective = %g)",
			        mpq_get_d (d_obj));
		}
		rval = QS_EXACT_UNSAT;
	}
	else
	{
		if(!msg_lvl)
		{
			MESSAGE(0, "Failed to prove problem unsatisfiable or delta-satisfiable (dual objective = %g)",
			        mpq_get_d (d_obj));
		}
		rval = QS_EXACT_UNKNOWN;
	}

	/* ending */
CLEANUP:
	mpq_EGlpNumFreeArray (dz);
	mpq_EGlpNumFreeArray (rhs_copy);
	mpq_clear (num1);
	mpq_clear (num2);
	mpq_clear (num3);
	mpq_clear (d_obj);
	return rval;
}

/* ========================================================================= */
/** @brief Used as separator while printing output to the screen (controled by
 * enabling simplex_display in the mpq_QSdata */
/* ========================================================================= */
static const char __sp[81] =
	"================================================================================";

/* ========================================================================= */
/** @brief print onto screen (if enable) a message indicating that we have
 * successfully solved the delta-satisfiability problem, and save (if x and y
 * are non NULL respectivelly) the optimal primal/dual solution provided in
 * x_mpq and y_mpq.
 * @param p_mpq the problem data.
 * @param x where to store the optimal primal solution (if not null).
 * @param y where to store the optimal dual solution (if not null).
 * @param x_mpq  the optimal primal solution.
 * @param y_mpq  the optimal dual solution.
 * @param delta  the tolerance parameter.
 * */
/* ========================================================================= */
static void delta_solved_output (mpq_QSdata * p_mpq,
																 mpq_t * const x,
																 mpq_t * const y,
																 mpq_t * x_mpq,
																 mpq_t * y_mpq,
																 int sat_status,
																 mpq_t const delta)
{
	if (p_mpq->simplex_display)
	{
		QSlog("delta-satisfiability problem solved exactly (result = %s)"
		      " with delta = %lg",
					sat_status == QS_EXACT_SAT ? "sat"
				: sat_status == QS_EXACT_UNSAT ? "unsat"
				: sat_status == QS_EXACT_DELTA_SAT ? "delta-sat" : "unknown",
				  mpq_get_d(delta));
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
/* Important: the input problem must be feasible and bounded.  */
int QSexact_delta_solver (mpq_QSdata * p_mpq,
													mpq_t * const x,
													mpq_t * const y,
													QSbasis * const ebasis,
													int simplexalgo,
													int *sat_status,
													mpq_t const delta)
{
	/* local variables */
	int status = 0, last_status = 0, last_iter = 0;
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
	*sat_status = QS_EXACT_UNKNOWN;
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
	EGcallD(dbl_QSget_status (p_dbl, &status));
	if ((status == QS_LP_INFEASIBLE) &&
			(p_dbl->lp->final_phase != PRIMAL_PHASEI) &&
			(p_dbl->lp->final_phase != DUAL_PHASEII))
		EGcallD(dbl_QSopt_primal (p_dbl, &status));
	EGcallD(dbl_QSget_status (p_dbl, &status));
	last_status = status;
	EGcallD(dbl_QSget_itcnt(p_dbl, 0, 0, 0, 0, &last_iter));
	/* optimization did not fail, so we have a (factorized) basis,
	 * which may be delta-sat */
	x_dbl = dbl_EGlpNumAllocArray (p_dbl->qslp->ncols);
	y_dbl = dbl_EGlpNumAllocArray (p_dbl->qslp->nrows);
	EGcallD(dbl_QSexact_delta_force_grab_cache (p_dbl, status, 1));
	EGcallD(dbl_QSget_x_array (p_dbl, x_dbl));
	EGcallD(dbl_QSget_pi_array (p_dbl, y_dbl));
	x_mpq = QScopy_array_dbl_mpq (x_dbl);
	y_mpq = QScopy_array_dbl_mpq (y_dbl);
	dbl_EGlpNumFreeArray (x_dbl);
	dbl_EGlpNumFreeArray (y_dbl);
	basis = dbl_QSget_basis (p_dbl);
	MESSAGE (msg_lvl, "Basis hash is 0x%016lX", QSexact_basis_hash(basis));
	*sat_status = QSexact_delta_optimal_test (p_mpq, x_mpq, y_mpq, basis, delta);
	if (QS_EXACT_UNKNOWN != *sat_status)
	{
		delta_solved_output (p_mpq, x, y, x_mpq, y_mpq, *sat_status, delta);
		goto CLEANUP;
	}
	MESSAGE (msg_lvl, "Retesting solution in exact arithmetic");
	EGcallD(QSexact_basis_status (p_mpq, &status, basis, msg_lvl, &simplexalgo));
	/* exact pivoting did not fail, so we have solution values,
	 * which may be delta-sat */
	EGcallD(mpq_QSexact_delta_force_grab_cache (p_mpq, status, 0));
	EGcallD(mpq_QSget_x_array (p_mpq, x_mpq));
	EGcallD(mpq_QSget_pi_array (p_mpq, y_mpq));
	*sat_status = QSexact_delta_optimal_test (p_mpq, x_mpq, y_mpq, basis, delta);
	if (QS_EXACT_UNKNOWN != *sat_status)
	{
		delta_solved_output (p_mpq, x, y, x_mpq, y_mpq, *sat_status, delta);
		goto CLEANUP;
	}
	else
	{
		last_status = status = QS_LP_UNSOLVED;
	}
	mpq_EGlpNumFreeArray (x_mpq);
	mpq_EGlpNumFreeArray (y_mpq);
	/* if we reach this point, then we have to keep going, we use the previous
	 * basis ONLY if the previous precision think that it has the optimal
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
		EGcallD(mpf_QSget_status (p_mpf, &status));
		if ((status == QS_LP_INFEASIBLE) &&
				(p_mpf->lp->final_phase != PRIMAL_PHASEI) &&
				(p_mpf->lp->final_phase != DUAL_PHASEII))
			EGcallD(mpf_QSopt_primal (p_mpf, &status));
		EGcallD(mpf_QSget_status (p_mpf, &status));
		last_status = status;
		EGcallD(mpf_QSget_itcnt(p_mpf, 0, 0, 0, 0, &last_iter));
		/* optimization did not fail, so we have a (factorized) basis,
		 * which may be delta-sat */
		basis = mpf_QSget_basis (p_mpf);
		MESSAGE(msg_lvl, "Basis hash is 0x%016lX", QSexact_basis_hash(basis));
		x_mpf = mpf_EGlpNumAllocArray (p_mpf->qslp->ncols);
		y_mpf = mpf_EGlpNumAllocArray (p_mpf->qslp->nrows);
		EGcallD(mpf_QSexact_delta_force_grab_cache (p_mpf, status, 1));
		EGcallD(mpf_QSget_x_array (p_mpf, x_mpf));
		EGcallD(mpf_QSget_pi_array (p_mpf, y_mpf));
		x_mpq = QScopy_array_mpf_mpq (x_mpf);
		y_mpq = QScopy_array_mpf_mpq (y_mpf);
		mpf_EGlpNumFreeArray (x_mpf);
		mpf_EGlpNumFreeArray (y_mpf);
		*sat_status = QSexact_delta_optimal_test (p_mpq, x_mpq, y_mpq, basis, delta);
		if (QS_EXACT_UNKNOWN != *sat_status)
		{
			delta_solved_output (p_mpq, x, y, x_mpq, y_mpq, *sat_status, delta);
			goto CLEANUP;
		}
		MESSAGE (msg_lvl, "Retesting solution in exact arithmetic");
		EGcallD(QSexact_basis_status (p_mpq, &status, basis, msg_lvl, &simplexalgo));
		/* exact pivoting did not fail, so we have solution values,
		 * which may be delta-sat */
		EGcallD(mpq_QSexact_delta_force_grab_cache (p_mpq, status, 0));
		EGcallD(mpq_QSget_x_array (p_mpq, x_mpq));
		EGcallD(mpq_QSget_pi_array (p_mpq, y_mpq));
		*sat_status = QSexact_delta_optimal_test (p_mpq, x_mpq, y_mpq, basis, delta);
		if (QS_EXACT_UNKNOWN != *sat_status)
		{
			delta_solved_output (p_mpq, x, y, x_mpq, y_mpq, *sat_status, delta);
			goto CLEANUP;
		}
		else
		{
			last_status = status = QS_LP_UNSOLVED;
		}
		mpq_EGlpNumFreeArray (x_mpq);
		mpq_EGlpNumFreeArray (y_mpq);
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

/* ========================================================================= */
/** @} */
/* end of exact_delta.c */

