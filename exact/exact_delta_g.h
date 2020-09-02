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
#ifndef EGLPNUM_TYPENAME__EXACT_DELTA_G_H__
#define EGLPNUM_TYPENAME__EXACT_DELTA_G_H__

#include "qstruct_EGLPNUM_TYPENAME.h"

/** @file
 * @ingroup Esolver */
/** @addtogroup Esolver */
/** @{ */
/* ========================================================================= */

/* ========================================================================= */
int EGLPNUM_TYPENAME_QSexact_delta_force_grab_cache (EGLPNUM_TYPENAME_QSdata *p,
																										 int status,
																										 int compute_pII_vals);

/** @} */
/* ========================================================================= */
/* end of exact_delta_g.h */
#endif

