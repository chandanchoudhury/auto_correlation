/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009,2010,2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include <gromacs/copyrite.h>
#include <gromacs/filenm.h>
#include <gromacs/macros.h>
#include <gromacs/pbc.h>
#include <gromacs/smalloc.h>
#include <gromacs/statutil.h>
#include <gromacs/xvgr.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/rmpbc.h>
#include <gromacs/nbsearch.h>
#include <gromacs/trajana.h>
#include <gromacs/typedefs.h>
#include <gromacs/vec.h>
#include <gromacs/tpxio.h>
#include <math.h>
#include <gromacs/index.h>
#include <gromacs/typedefs.h>
#include<time.h>
#include"varS.h"
/*! \brief
 * Template analysis data structure.
 */

int main(int argc, char *argv[]) {
   clock_t start = clock();
   const char *desc[] = {
       "this is a small test program meant to serve as a template ", 
	   "when writing your own analysis tools. The advantage of this",
	   "is you can have access to the gromcas parameters in the",
	   "topology and also you can read the binary trajectories" 
	};
	
	int n = 1;
	static gmx_bool bPBC         = FALSE;
	t_pargs pa[] = {
    	{ "-n", FALSE, etINT, {&n},
      	"Plot data for atom number n (starting on 1)" },
		{ "-pbc",      FALSE, etBOOL, {&bPBC},
		"Use periodic boundary conditions for computing distances" }
	};
	
	t_topology 		top;
	int        		ePBC;
	char       		title[STRLEN];
	t_trxframe 		fr;
	rvec       		*xtop;
	matrix     		box, box_pbc;
	int        		flags = TRX_READ_X;
	t_trxstatus 	*status;
	t_pbc       	pbc;
    t_atoms    		*atoms;
    int         	natoms;
	char       		*grpnm, *grpnmp;
    atom_id    		*index, *indexp;
    int         	i, nidx, nidxp;
    int        		v;
    int         	j, k;
    long         	***bin = (long ***)NULL;
	output_env_t    oenv;
	gmx_rmpbc_t     gpbc=NULL;

	t_filenm fnm[] = {
	    { efTPS,  NULL,  NULL, ffREAD },   /* this is for the topology */
    	{ efTRX, "-f", NULL, ffREAD }      /* and this for the trajectory */
  	};

	#define NFILE asize(fnm)

	parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);


	read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
	sfree(xtop);

	n=n-1; /* Our enumeration started on 1, but C starts from 0 */
  	/* check that this atom exists */
  	if(n<0 || n>(top.atoms.nr)) {
	    printf("Error: Atom number %d is out of range.\n",n);
    	exit(1);
  	}
	
 /* RESINFO are meant for pdb;
 	printf("Atom name: %s\n",*(top.atoms.atomname[n]));
  	printf("Atom charge: %f\n",top.atoms.atom[n].q);
  	printf("Atom Mass: %f\n",top.atoms.atom[n].m);
  	printf("Atom resindex: %d\n",top.atoms.atom[n].resind); //Starts from 0
  	printf("Atom resname: %s\n",*(top.atoms.resinfo[top.atoms.atom[n].resind].name));
  	printf("Atom resno: %d %d\n",i, (top.atoms.resinfo[top.atoms.atom[i].resind].nr));
  	printf("Atom chain id: %d %c\n",i, top.atoms.resinfo[top.atoms.atom[i].resind].chainid); //Represents the pdb numbering
  	printf("Atom chain no: %d %d\n",i,top.atoms.resinfo[top.atoms.atom[i].resind].chainnum);
  	printf("Atom resno: %d %d\n",i, top.atoms.atom[i].resind);
*/
	atoms = &(top.atoms);

	/* The first time we read data is a little special */
	natoms = read_first_frame(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &fr, flags);

	/*My Code Starts Here*/

	FILE *facf_dis;
	facf_dis = fopen("acfS.xvg","wb");
    int is = 0, ie = 0, chain = 0, d = 0, frame = 0;
    int row = 60000, col = tot_chains + 2, f = 0;
    double distsq = 0.0, **R, Rsumall = 0.0, Ravg = 0.0;
    double R2sumall = 0.0, normalize = 1.0, std = 0.0, min = +9.999e99, max = -9.999e99;
    double R2sum = 0.0, acf[coor_length];

    memset(acf,   0, coor_length * 1 * sizeof(double));

    R = (double**) malloc ( sizeof(double*) * row);
    if (R == NULL) {
        printf("*** ERROR - Could not malloc() enough space \n"); exit(1);
    }
    memset(R, 0, sizeof(double) * (row));

    for(i = 0; i < row; i++) {
        R[i] = (double *) malloc ( sizeof(double) * col);
        memset(R[i], 0, sizeof(double) * col);
    }

    fprintf(stdout,"\n\n");
    fprintf(stdout,"Coorelation Length           = %4d\n", coor_length);
    fprintf(stdout,"Total Chains in System       = %4d\n", tot_chains);
    fprintf(stdout,"Beads per Chain              = %4d\n", beads_per_chain);
    fprintf(stdout,"Frame to Time(ps) Conversion = %4d\n", frame2ps);
    fprintf(stdout,"\n\n");

    do {
        is = 0; ie = 0; chain = 0; distsq = 0.0;

        for(chain = 1; chain <= tot_chains; chain++) {
            ie = chain * beads_per_chain - 1;
            distsq = 0.0;
        
            for(d = 0; d < DIM; d++) 
                distsq += pow((fr.x[is][d] - fr.x[ie][d]),2);

            R[frame][chain] = sqrt(distsq);
            Rsumall += R[frame][chain];

            R2sumall += distsq;

//            max = (max > R[frame][chain]) ? max : R[frame][chain];
//            min = (min < R[frame][chain]) ? min : R[frame][chain];

            is = ie + 1;
        }
        frame++;
    }
    while (read_next_frame(oenv, status, &fr));

    normalize = 1.0/(double)(frame * tot_chains);

    Ravg = Rsumall * normalize;
    std  = sqrt((R2sumall * normalize) - (Ravg * Ravg));

//    fprintf(stderr,"   Avg        Std        Max         Min        Frames\n");
//    fprintf(stderr,"%10g %10g %10g %10g %10d\n", Ravg,  std, max, min, frame);

    std = 1.0/std;

    for(f = 0; f < frame; f++)
        for(chain = 1; chain <= tot_chains; chain++)
            R[f][chain] = (R[f][chain] - Ravg) * std;
            
    for(i = 0; i < coor_length; i++) {
        acf[i] = 0.0;
        for (f = 0; f < (frame - i); f++) 
            for(chain = 1; chain <= tot_chains; chain++) 
                acf[i] += R[f][chain] * R[f + i][chain];

        acf[i] /= ((double) (tot_chains * (frame - i)));
        fprintf(facf_dis,"%4d %8.4f\n",i * frame2ps, acf[i]);
    }

    for(i = 0; i < row; i++)
        free(R[i]);
    
    free(R);
    printf("\nTime elapsed: %f s\n", ((double)clock() - start) / CLOCKS_PER_SEC);
}
