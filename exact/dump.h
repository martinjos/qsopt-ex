/* dump - Functions for outputting internal data structures
 *
 * Copyright (C) 2020  Martin Sidaway
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * */

#ifndef QS_DUMP_H__
#define QS_DUMP_H__

#include "qstruct_mpq.h"

void mpq_QSdump_xbz (const mpq_QSdata *p_mpq);
void mpq_QSdump_xnbz (const mpq_QSdata *p_mpq);
void mpq_QSdump_bfeas (const mpq_QSdata *p_mpq);
void mpq_QSdump_array (const mpq_t *array, const char* tag);

void mpq_QSdump_prob_col (const mpq_QSdata *p_mpq, int index, int col, char type);
void mpq_QSdump_prob (const mpq_QSdata *p_mpq);

// p_mpq can't be const because it is passed to mpq_ILLlib_tableau()
int mpq_QSdump_basis (mpq_QSdata *p_mpq);

#endif /* ! QS_DUMP_H__ */
