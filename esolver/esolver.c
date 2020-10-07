/* ========================================================================= */
/* ESolver "Exact Mixed Integer Linear Solver" provides some basic structures 
 * and algorithms commons in solving MIP's
 *
 * Copyright (C) 2008 David Applegate, Bill Cook, Sanjeeb Dash, Daniel Espinoza.
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
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <signal.h>
#include <limits.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "QSopt_ex.h"

#include "except.h"
#include "logging-private.h"
#include "qs_config.h"

enum {
	ALGO_OPT = 0,
	ALGO_DELTA_OPT,
	ALGO_DELTA_FEAS_PHASE_1,
	ALGO_DELTA_FEAS_PHASE_2,
};

/* ========================================================================= */
/** @name static parameters for the main program */
/*@{*/
static char *fname = 0;
static int lpfile = 0;
static int algo = ALGO_OPT;
static int usescaling = 1;
static int showversion = 0;
static int simplex_display = 1;
static int simplexalgo = PRIMAL_SIMPLEX;
static int pstrategy = QS_PRICE_PSTEEP;
static int dstrategy = QS_PRICE_DSTEEP;
static unsigned precision = 128;
static int printsol = 0;
static int continuous_output = 0;
static char *solname = 0;
static char *readbasis = 0;
static char *writebasis = 0;
/** @brief maximum running time */
static double max_rtime = INT_MAX;
/** @brief maximum memory usage */
static unsigned long memlimit = UINT_MAX;
static mpq_t delta;
/*@}*/
/* ========================================================================= */
/** @brief Display options to the screen */
static void usage (char *s)
{
	fprintf (stderr, "Usage: %s [- below -] prob_file\n", s);
	fprintf (stderr, "   -a o,do,dfp1,dfp2  select algorithm\n");
	fprintf (stderr, "       o = optimality [default]\n");
	fprintf (stderr, "       do = delta-optimality\n");
	fprintf (stderr, "       dfp1 = delta-feasibility, phase 1\n");
	fprintf (stderr, "       dfp2 = delta-feasibility, phase 2\n");
	fprintf (stderr, "   -D #  delta to use for delta-complete algorithms (default: 0)\n");
	fprintf (stderr, "   -c    display all delta-sat thresholds when available (-a dfp1 or dfp2 only)\n");
	fprintf (stderr, "   -b f  write basis to file f\n");
	fprintf (stderr, "   -B f  read initial basis from file f\n");
#if 0
	fprintf (stderr, "   -I    solve the MIP using BestBound\n");
	fprintf (stderr, "   -E    edit problem after solving initial version\n");
#endif
	fprintf (stderr, "   -L    input file is in lp format (default: mps)\n");
	fprintf (stderr, "   -O    write the final solution to the given file\n");
	fprintf (stderr, "         append .gz/.bz2 to the .sol extension to compress the file\n");
	fprintf (stderr, "   -M    write model to stdout\n");
	fprintf (stderr, "   -p #  run primal simplex with pricing rule #\n");
	fprintf (stderr,
					 "         (%d-Dantzig, %d-Devex, %d-Steep (default), %d-Partial\n",
					 QS_PRICE_PDANTZIG, QS_PRICE_PDEVEX, QS_PRICE_PSTEEP,
					 QS_PRICE_PMULTPARTIAL);
	fprintf (stderr,
					 "   -P #  number of bits to use for the float representation (default: 128)\n");
	fprintf (stderr, "   -d #  run dual simplex with pricing rule #\n");
	fprintf (stderr, "         (%d-Dantzig, %d-Steep, %d-Partial, %d-Devex)\n",
					 QS_PRICE_DDANTZIG, QS_PRICE_DSTEEP, QS_PRICE_DMULTPARTIAL,
					 QS_PRICE_DDEVEX);
	fprintf (stderr, "   -S    do NOT scale the initial LP\n");
	fprintf (stderr, "   -v    print QSopt version number\n");
	fprintf (stderr, "   -V n  verbosity level (0-3, default: 1)\n");
	fprintf (stderr, "   -R n  maximum running time allowed, default %lf\n",
						max_rtime);
	fprintf (stderr, "   -m n  maximum memory usage allowed, default %lu\n", 
						memlimit);
}

/* ========================================================================= */
/** @brief decide if a given file is mps or lp (only by extension) */
static void get_ftype (char const *const name,
											 int *ftype)
{
	char buff[4096],*argv[128];
	int argc;
	*ftype = 0; /* by default, file is MPS */
	snprintf(buff,4096,"%s",name);
	EGioNParse(buff,128,"."," ",&argc,argv);
	argc-=1;
	if(argc)
	{
		if(strncmp(argv[argc],"gz",3)==0) argc-=1;
		else if(strncmp(argv[argc],"GZ",3)==0) argc-=1;
		else if(strncmp(argv[argc],"bz2",4)==0) argc-=1;
		else if(strncmp(argv[argc],"BZ2",4)==0) argc-=1;
	}
	if(argc)
	{
		if(strncmp(argv[argc],"lp",3)==0) *ftype=1;
		else if(strncmp(argv[argc],"LP",3)==0) *ftype=1;
	}
}
/* ========================================================================= */
/** @brief signal handler for time-limit reached */
static void sighandler(int s)
{
	switch(s)
	{
		case SIGXCPU:
			fprintf(stderr,"TIME_LIMIT_REACHED (ending now)\n");
			exit(EXIT_FAILURE);
		default:
			fprintf(stderr,"Unknown signal %d (ending now)\n",s);
			exit(EXIT_FAILURE);
	}
}
/* ========================================================================= */
/** @brief function to handle resource usage limits */
static int mem_limits(int log_level)
{
	int rval = 0;
	struct rlimit mlim;
	rval = getrlimit(RLIMIT_CPU,&mlim);
	CHECKRVAL(rval);
	if (log_level >= 1)
		fprintf(stderr, "Cur rtime limit %ld, trying to set to %lg\n", mlim.rlim_cur, max_rtime);
	if(max_rtime > mlim.rlim_max) max_rtime = (double)mlim.rlim_max;
	mlim.rlim_cur = (rlim_t)max_rtime;
	rval = setrlimit(RLIMIT_CPU,&mlim);
	TESTERRNOIF(rval);
	if (log_level >= 1)
		fprintf(stderr, "New rtime limit %ld (%.3lg)\n", mlim.rlim_cur, max_rtime);
	rval = getrlimit(RLIMIT_DATA,&mlim);
	TESTERRNOIF(rval);
	if (log_level >= 1)
		fprintf(stderr, "Cur data limit %ld,%ld (soft,hard)\n", mlim.rlim_cur,
						mlim.rlim_max);
	mlim.rlim_cur = memlimit;				
	rval = setrlimit(RLIMIT_DATA,&mlim);				
	TESTERRNOIF(rval);
	rval = getrlimit(RLIMIT_DATA,&mlim);
	TESTERRNOIF(rval);
	if (log_level >= 1)
		fprintf(stderr, "New data limit %ld,%ld (soft,hard)\n", mlim.rlim_cur,
						mlim.rlim_max);
	rval = getrlimit(RLIMIT_AS,&mlim);
	TESTERRNOIF(rval);
	if (log_level >= 1)
		fprintf(stderr, "Cur address space limit %ld,%ld (soft,hard)\n",
						mlim.rlim_cur, mlim.rlim_max);
	mlim.rlim_cur = memlimit;	
	rval = setrlimit(RLIMIT_AS,&mlim);				
	TESTERRNOIF(rval);
	rval = getrlimit(RLIMIT_AS,&mlim);
	TESTERRNOIF(rval);
	if (log_level >= 1)
		fprintf(stderr, "New address space limit %ld,%ld (soft,hard)\n",
						mlim.rlim_cur, mlim.rlim_max);
	mlim.rlim_cur = 0;
	rval = setrlimit(RLIMIT_CORE,&mlim);				
	TESTERRNOIF(rval);
	rval = getrlimit(RLIMIT_CORE,&mlim);
	TESTERRNOIF(rval);
	if (log_level >= 1)
		fprintf(stderr, "New core dump space limit %ld,%ld (soft,hard)\n",
						mlim.rlim_cur, mlim.rlim_max);
	/* set signal handler for SIGXCPU */
	signal(SIGXCPU,sighandler);
	return rval;
}
/* ========================================================================= */
/** @brief parssing options for the program */
static int parseargs (int ac,
											char **av)
{
	int c;
	int boptind = 1;
	char *boptarg = 0;

	while ((c = ILLutil_bix_getopt (ac, av, "a:b:cB:d:D:EILm:MO:p:P:R:SvV:",
																	&boptind, &boptarg)) != EOF)
		switch (c)
		{
		case 'a':
			if (!strcmp("o", boptarg)) {
				algo = ALGO_OPT;
			} else if (!strcmp("do", boptarg)) {
				algo = ALGO_DELTA_OPT;
			} else if (!strcmp("dfp1", boptarg)) {
				algo = ALGO_DELTA_FEAS_PHASE_1;
			} else if (!strcmp("dfp2", boptarg)) {
				algo = ALGO_DELTA_FEAS_PHASE_2;
			} else {
				usage (av[0]);
				return 1;
			}
			break;
		case 'c':
			continuous_output = 1;
			break;
		case 'D': {
			int nchar = mpq_EGlpNumReadStrXc (delta, boptarg);
			if (nchar != strlen(boptarg)) {
				fprintf (stderr, "Failed to parse number: %s\n\n", boptarg);
				usage (av[0]);
				return 1;
			}
			break;
		}
		case 'm':
			memlimit = strtoul(boptarg,0,10);
			break;
		case 'R':
			max_rtime = strtod(boptarg,0);
			break;
		case 'b':
			writebasis = boptarg;
			break;
		case 'B':
			readbasis = boptarg;
			break;
		case 'P':
			precision = atoi (boptarg);
			break;
		case 'd':
			simplexalgo = DUAL_SIMPLEX;
			dstrategy = atoi (boptarg);
			break;
		case 'L':
			lpfile = 1;
			break;
		case 'O':
			printsol = 1;
			solname = strdup(boptarg);
			break;
		case 'M':
			printsol = 1;
			break;
		case 'p':
			simplexalgo = PRIMAL_SIMPLEX;
			pstrategy = atoi (boptarg);
			break;
		case 'S':
			usescaling = 0;
			break;
		case 'v':
			showversion = 1;
			break;
		case 'V':
			simplex_display = atoi (boptarg);
			break;
		case '?':
		default:
			usage (av[0]);
			return 1;
		}
	if ((boptind == ac) && (showversion))
	{
		char *buf = 0;
		buf = mpq_QSversion ();
		printf ("%s\n", buf);
		mpq_QSfree ((void *) buf);
		exit(0);
	}
	if (boptind != (ac - 1))
	{
		usage (av[0]);
		return 1;
	}

	fname = av[boptind++];
	fprintf (stderr, "Reading problem from %s\n", fname);
	mem_limits(simplex_display);
	return 0;
}

void free_gmp_str (char *str)
{
	// Free memory allocated by mpq_get_str (etc.)
	void (*free_func) (void *, size_t);
	mp_get_memory_functions (NULL, NULL, &free_func);
	free_func (str, strlen (str) + 1);
}

void partial_result (mpq_QSdata const * p_mpq,
										 mpq_t * const x,
										 const mpq_t infeas,
										 const mpq_t delta,
										 void *data)
{
	char *mpq_str = NULL;
	mpq_str = mpq_get_str (NULL, 10, infeas);
	// Strictly speaking, this is > rather than >=, but the point is that it
	// could be anything >=.
	printf ("PARTIAL: status = delta-FEASIBLE with delta = %.17g ( >= %s)\n",
					nextafter (mpq_get_d (infeas), 1.0/0.0), mpq_str);
	free_gmp_str (mpq_str);
}

/* ========================================================================= */
/** @brief the main thing! */
/* ========================================================================= */
int main (int ac,
					char **av)
{
	int rval = 0,
	  status = 0;
	mpq_QSdata *p_mpq = 0;
	QSbasis *basis = 0;
	ILLutil_timer timer_solve;
	ILLutil_timer timer_read;
	int ftype = 0;								/* 0 mps, 1 lp */
	mpq_t *y_mpq = 0,
	 *x_mpq = 0;
	mpq_t obj_lo;
	mpq_t obj_up;
	QSopt_ex_version();
	QSexactStart();
	mpq_init(delta);  // init to 0
	mpq_init(obj_lo);
	mpq_init(obj_up);
	/* parse arguments and initialize EGlpNum related things */
	rval = parseargs (ac, av);
	QSexact_set_precision (precision);
	if (rval)
		goto CLEANUP;
	if (writebasis)
	{
		basis = EGsMalloc (QSbasis, 1);
		memset (basis, 0, sizeof (QSbasis));
	}

	/* just for the bell's and wistle */
	if (showversion)
	{
		char *buf = 0;
		buf = mpq_QSversion ();
		if (buf == 0)
		{
			ILL_CLEANUP;
		}
		else
		{
			printf ("%s\n", buf);
			mpq_QSfree ((void *) buf);
		}
	}

	/* get the file type */
	if (lpfile)
		ftype = 1;
	else
		get_ftype (fname, &ftype);

	/* read the mpq problem */
	ILLutil_init_timer (&timer_read, "SOLVER_READ_MPQ");
	ILLutil_start_timer (&timer_read);
	if (ftype == 1)
	{
		p_mpq = mpq_QSread_prob ((const char *) fname, "LP");
		if (p_mpq == 0)
		{
			fprintf (stderr, "Could not read lp file.\n");
			rval = 1;
			ILL_CLEANUP_IF (rval);
		}
	}
	else
	{
		p_mpq = mpq_QSread_prob ((const char *) fname, "MPS");
		if (p_mpq == 0)
		{
			fprintf (stderr, "Could not read mps file.\n");
			rval = 1;
			ILL_CLEANUP_IF (rval);
		}
	}

	/* and get the basis if needed */
	if (readbasis)
	{
		rval = mpq_QSread_and_load_basis (p_mpq, (const char *) readbasis);
		ILL_CLEANUP_IF (rval);
		if (basis)
			mpq_QSfree_basis (basis);
		basis = mpq_QSget_basis (p_mpq);
	}
	ILLutil_stop_timer (&timer_read, 1);
	/* set the readed flags */
	rval = mpq_QSset_param (p_mpq, QS_PARAM_SIMPLEX_DISPLAY, simplex_display)
		|| mpq_QSset_param (p_mpq, QS_PARAM_PRIMAL_PRICING, pstrategy)
		|| mpq_QSset_param (p_mpq, QS_PARAM_DUAL_PRICING, dstrategy)
		|| mpq_QSset_param (p_mpq, QS_PARAM_SIMPLEX_SCALING, usescaling);
	ILL_CLEANUP_IF (rval);
	if (printsol || algo == ALGO_DELTA_OPT)
	{
		x_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->ncols);
		y_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->nrows);
	}
	ILLutil_init_timer (&timer_solve, "SOLVER");
	ILLutil_start_timer (&timer_solve);
	switch (algo)
	{
	case ALGO_OPT:
		rval = QSexact_solver (p_mpq, x_mpq, y_mpq, basis, simplexalgo, &status);
		break;
	case ALGO_DELTA_OPT:
		rval = QSdelta_full_solver (p_mpq, delta, x_mpq, y_mpq, obj_lo, obj_up, basis, simplexalgo, &status);
		mpq_sub (delta, obj_up, obj_lo);
		break;
	case ALGO_DELTA_FEAS_PHASE_1:
		rval = QSdelta_solver (p_mpq, delta, x_mpq, y_mpq, basis, simplexalgo, &status,
													 continuous_output ? partial_result : NULL, NULL);
		break;
	case ALGO_DELTA_FEAS_PHASE_2:
		rval = QSexact_delta_solver (p_mpq, x_mpq, y_mpq, basis, simplexalgo, &status, delta,
																 continuous_output ? partial_result : NULL, NULL);
		break;
	}
	ILL_CLEANUP_IF (rval);
	ILLutil_stop_timer (&timer_solve, 1);
	/* display result */
	EGioFile_t *out_f;
	if (solname)
	{
		char out_f_name[1024];
		sprintf (out_f_name, "%s", solname);
		out_f = EGioOpen (out_f_name, "w");
	}
	else
		out_f = EGioOpenFILE (stdout);
	switch (status)
	{
	case QS_LP_OPTIMAL:
		EGioPrintf (out_f, "status = OPTIMAL\n");
		if (printsol)
		{
			rval = QSexact_print_sol (p_mpq, out_f);
			CHECKRVALG(rval,CLEANUP);
		}
		break;
	case QS_LP_DELTA_OPTIMAL:
		EGioPrintf (out_f, "status = delta-OPTIMAL with delta = %g\n", delta);
		if (printsol)
		{
			rval = QSexact_print_sol (p_mpq, out_f);
			CHECKRVALG(rval,CLEANUP);
		}
		break;
	case QS_LP_FEASIBLE:
		EGioPrintf (out_f, "status = FEASIBLE\n");
		if (printsol)
		{
			rval = QSexact_print_sol (p_mpq, out_f);
			CHECKRVALG(rval,CLEANUP);
		}
		break;
	case QS_LP_DELTA_FEASIBLE:
		EGioPrintf (out_f, "status = delta-FEASIBLE with delta = %g\n", delta);
		if (printsol)
		{
			rval = QSexact_print_sol (p_mpq, out_f);
			CHECKRVALG(rval,CLEANUP);
		}
		break;
	case QS_LP_INFEASIBLE:
		EGioPrintf (out_f, "status = INFEASIBLE\n");
		break;
	case QS_LP_UNBOUNDED:
		EGioPrintf (out_f, "status = UNBOUNDED\n");
		break;
	default:
		EGioPrintf (out_f, "status = UNDEFINED\n");
		break;
	}
	EGioClose (out_f);
	/* ending */
CLEANUP:
	if (solname)
		EGfree(solname);
	mpq_EGlpNumFreeArray (x_mpq);
	mpq_EGlpNumFreeArray (y_mpq);
	/* free the last allocated basis, and if we wanted to save it, do so */
	if (basis)
	{
		if (writebasis)
			rval = mpq_QSwrite_basis (p_mpq, 0, writebasis);
	}
	mpq_QSfree_basis (basis);
	mpq_QSfree_prob (p_mpq);
	mpq_clear(obj_up);
	mpq_clear(obj_lo);
	mpq_clear(delta);
	QSexactClear();
	return rval;									/* main return */
}
