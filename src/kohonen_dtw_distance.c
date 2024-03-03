//
// Copyright (c) 2006-2019 of Toni Giorgino
//
// This file is part of the DTW package.
//
// DTW is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DTW is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with DTW.  If not, see <http://www.gnu.org/licenses/>.
//

// If you copy this algorithm, please cite doi:10.18637/jss.v031.i07

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "./kohonen_dtw_distance.h"

// Define R-like functions - Memory alloc'd by R_alloc is automatically freed
#ifdef DTW_R
#include <R.h>
#define dtw_alloc(n, size) (R_alloc(n, size))
#else
// Either standalone or Python
#include <limits.h>
#define R_NaInt INT_MIN
#if _WIN32
#include <malloc.h>
#define dtw_alloc(n, size) (_alloca((n) * (size)))
#else
#define dtw_alloc(n, size) (alloca((n) * (size)))
#endif
#endif

#ifndef NAN
#error "This code requires native IEEE NAN support. Verify you are using gcc with -std=gnu99 or recent compilers."
#endif

/* undo R indexing */
#define EP(ii, jj) ((jj) * nsteps + (ii))
#define EM(ii, jj) ((jj) * n + (ii))

#define CLEARCLIST                      \
    {                                   \
        for (int z = 0; z < npats; z++) \
            clist[z] = NAN;             \
    }

/*
 * Auxiliary function: return the arg min, ignoring NANs, -1 if all NANs
 * TODO: remove isnan and explain, check, time
 */
static inline int argmin(const double *list, int n)
{
    int ii = -1;
    double vv = INFINITY;
    for (int i = 0; i < n; i++)
    {
        /* The following is a faster equivalent to
         *    if(!isnan(list[i]) && list[i]<vv)
         * because   (NAN < x) is false
         */
        if (list[i] < vv)
        {
            ii = i;
            vv = list[i];
        }
    }
    return ii;
}

/*
 *  Compute cumulative cost matrix: replaces kernel in globalCostMatrix.R
 */

/* For now, this code is also valid outside R, as a test unit.
 * This means that we have to refrain to use R-specific functions,
 * such as R_malloc, or conditionally provide replacements.
 */

/* R matrix fastest index is row */

void computeCM(        /* IN */
   const int s,        /* mtrx dimensions, int */
   const int *wm,      /* windowing matrix, logical=int */
   const double *lm,   /* local cost mtrx, numeric */
   const int *nstepsp, /* no of steps in stepPattern, int */
   const double *dir,  /* stepPattern description, numeric */
   /* IN+OUT */
   double *cm /* cost matrix, numeric */
   /* OUT - sits: Disabled as we do not use it */
   // int *sm /* direction mtrx, int */
)
{
    /* recover matrix dim */
    int n = s, m = s; /* query,template as usual*/
    int nsteps = *nstepsp;

    /* copy steppattern description to ints,
     so we'll do indexing arithmetic on ints
     */
    int *pn, *di, *dj;
    double *sc;

    pn = (int *)dtw_alloc((size_t)nsteps, sizeof(int));       /* pattern id */
    di = (int *)dtw_alloc((size_t)nsteps, sizeof(int));       /* delta i */
    dj = (int *)dtw_alloc((size_t)nsteps, sizeof(int));       /* delta j */
    sc = (double *)dtw_alloc((size_t)nsteps, sizeof(double)); /* step cost */

    for (int i = 0; i < nsteps; i++)
    {
        pn[i] = (int)dir[EP(i, 0)] - 1; /* Indexing C-way */
        di[i] = (int)dir[EP(i, 1)];
        dj[i] = (int)dir[EP(i, 2)];
        sc[i] = dir[EP(i, 3)];
    }

    /* assuming pattern ids are in ascending order */
    int npats = pn[nsteps - 1] + 1;

    /* prepare a cost list per pattern */
    double *clist = (double *)
        dtw_alloc((size_t)npats, sizeof(double));

    /* we do not initialize the seed - the caller is supposed
     to do so
     cm[0]=lm[0];
     */

    // sits: Disabled as we do not use it
    /* clear the direction matrix */
    // for (int i = 0; i < m * n; i++)
    // sm[i] = R_NaInt; /* should be NA_INTEGER? */

    /* lets go */
    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {

            /* out of window? */
            if (!wm[EM(i, j)])
                continue;

            /* already initialized? */
            if (!isnan(cm[EM(i, j)]))
                continue;

            CLEARCLIST;
            for (int s = 0; s < nsteps; s++)
            {
                int p = pn[s]; /* indexing C-way */

                int ii = i - di[s];
                int jj = j - dj[s];

                if (ii >= 0 && jj >= 0)
                { /* address ok? C convention */
                    double cc = sc[s];
                    if (cc == -1.0)
                    {
                        clist[p] = cm[EM(ii, jj)];
                    }
                    else
                    { /* we rely on NAN to propagate */
                        clist[p] += cc * lm[EM(ii, jj)];
                    }
                }
            }

            int minc = argmin(clist, npats);
            if (minc > -1)
            {
                cm[EM(i, j)] = clist[minc];
                // sits: Disabled as we do not use it
                // sm[EM(i, j)] = minc + 1; /* convert to 1-based  */
            }
        }
    }
}
