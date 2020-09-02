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

#include "exact_delta_g_EGLPNUM_TYPENAME.h"

#include "logging-private.h"

#include "except.h"

#include "lpdata_EGLPNUM_TYPENAME.h"
#include "fct_EGLPNUM_TYPENAME.h"

int EGLPNUM_TYPENAME_QSexact_delta_force_grab_cache (EGLPNUM_TYPENAME_QSdata *p,
																										 int status,
																										 int compute_pII_vals)
{
	int rval = 0;
	char save_optimal = p->lp->basisstat.optimal;
	if (status != QS_LP_OPTIMAL)
	{
		if (compute_pII_vals)
		{
			// Only if not already done by QSexact_basis_status
			EGLPNUM_TYPENAME_ILLfct_compute_xbz (p->lp);
			EGLPNUM_TYPENAME_ILLfct_compute_piz (p->lp);
			EGLPNUM_TYPENAME_ILLfct_compute_dz (p->lp);
		}
		p->lp->basisstat.optimal = 1;  // Pretend that it is optimal
		EGcallD(EGLPNUM_TYPENAME_QSgrab_cache (p, status));
	}
CLEANUP:
	p->lp->basisstat.optimal = save_optimal;
	EG_RETURN (rval);
}

/* ========================================================================= */
/** @} */
/* end of exact_delta_g.c */

