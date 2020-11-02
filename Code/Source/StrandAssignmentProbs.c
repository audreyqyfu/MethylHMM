/*****************************************************************************************
 * StrandAssignmentProbs.c
 *
 * 2010.09.09: Discrete Markov model representation.
 * 2010.04.16: allow rates of maintenance and de novo activities of DNMT3 to be different
 * 2009.08.09: parametrization as in Note_10_EmissProbs
 * July 26, 2009: file name and function name changed.
 * September 7, 2008: compute and return strand type probs
 * August 18, 2008: revised emission probabilities (see hmm_July2808)
 *
 * Function to compute the probability of assigning the top strand 
 * in a double-stranded sequence to be the parent strand under the 
 * hidden Markov model for processivity.  The Markov processces are 
 * modelled as jump processes with bisulfite conversion rates b 
 * (failure of conversion) and c (inappropriate conversion) constant 
 * across sites.
 * 
 * lambda contains 3 pi's, Pr(0->1) for M, RP, RD, vector of length 3*(NSite-1); 
 * Pr(1->0) for M, RP, RD, vector of length 3*(NSite-1); m, vector of length NSite
 *
 * input: data: unprocessed double-stranded data
 *        lambda: vector of length 6+3*2*(NSite-1)+NSite
 *		  MDenovo: de novo rate of Dnmt1
 *		  MMaint: rate of maintenance activity of Dnmt1
 *		  DDenovo: de novo rate of Dnmt3 on daughter strand
 *		  DMaint: rate of maintenance activity of Dnmt3 on daughter strand
 *        NSeq: number of ds seqs
 *        NSite: number of CpG sites
 *        NHidden: number of hidden states
 *        NObs: number of observed states
 *        DE_NOVO: 0 for assumption 1; 1 for assumption 2
 *        PRINT_FLAG: 0 for no output; 1 for output
 * output: loglik: log likelihood of all the data
 *		   strandtype: probs of top strand being parent; vector of NSite
 ******************************************************************************************/
# include <R.h>
# include <stdlib.h>
# include <stdio.h>
# include <math.h>

void compute_pi (double * pi, double * pi_ind);  
void compute_trans_mtx_jump (double *** trans_mtx, double ** rho_0, double ** rho_1, int n_proc, int n_site, int print_flag);
void compute_error_mtx (double ** error_mtx, double b, double c, int print_flag);
void compute_emiss_mtx_c (double *** emiss_mtx_c, double * m, double m_denovo, double m_maint, double d_denovo, double d_maint, double ** error_mtx, int n_site, int n_hidden, int n_obs, int DE_NOVO, int print_flag);
void multiply_mtx (double ** mtx_prod, double ** mtx_1, double ** mtx_2, int dim_i, int dim_j, int dim_k, int print_flag);
void Kronecker_prod (double ** mtx1, int dim1, double ** mtx2, int dim2, double ** mtx);
double loglikind_strandtype_hmm (int ** ind_seq, double * pi, double *** A, double *** B, int length, int n_hidden, int print_flag);
double likind_ordered_hmm (int ** ind_seq, double * pi, double *** A, double *** B, int length, int nHiddenState, int print_flag);
double forward_backward (double * pi, double *** A, double *** B, int * Obs, int length, int nHiddenState, int print_flag);


void StrandAssignmentProbs (double * strandtype, int * data, double * lambda, double * MDenovo, double * MMaint, double * DDenovo, double * DMaint, double * b, double * c, int * NSeq, int * NSite, int * NHidden, int * NObs, int * NProc, int * DE_NOVO, int * PRINT_FLAG)
{
  double m_denovo = MDenovo[0];		/* de novo rate of DNMT1 */
  double m_maint = MMaint[0];			/* rate of maintenance activity of DNMT1 */
	double d_denovo = DDenovo[0];		/* de novo rate of DNMT3s on daughter strand */
	double d_maint = DMaint[0];			/* rate of maintenance activity of DNMT3s on daughter strand */
  double error_b = b[0];
  double error_c = c[0];
  int de_novo_flag = DE_NOVO[0];     /* 0 for assumption 1; 1 for assumption 2 */
  int print_flag = PRINT_FLAG[0];    /* 0 for no output; 1 for output */
  int n_seq = NSeq[0];
  int n_site = NSite[0];
  int n_hidden = NHidden[0];
  int n_obs = NObs[0];
  int n_proc = NProc[0];         /* number of meth. processes; 3 here */
	double result = 0.0;
	double ** prob_0to1, ** prob_1to0;
  double * pi_ind, * pi, *** trans_mtx, *** emiss_mtx;
  double ** error_mtx;
  int ** ind_seq;
  double * m;
  double m_s;
  int i,j,k,s;

  /*************************************
   * allocate memory
   ************************************/
  pi = (double *) malloc (sizeof (double) * n_hidden);
  trans_mtx = (double ***) malloc (sizeof (double **) * (n_site-1));
  for (s=0; s<n_site-1; s++)
    {
      trans_mtx[s] = (double **) malloc (sizeof(double*)*n_hidden);
      for (i=0; i<n_hidden; i++)
		  trans_mtx[s][i] = (double *) malloc(sizeof(double)*n_hidden);
    }
  emiss_mtx = (double ***) malloc(sizeof(double**)*n_site);
  for (s=0; s<n_site; s++)
    {
      emiss_mtx[s] = (double **) malloc (sizeof(double*)*n_hidden);
      for (i=0; i<n_hidden; i++)
	emiss_mtx[s][i] = (double *) malloc(sizeof(double)*n_obs);
    }
  error_mtx = (double **) malloc (sizeof (double *) * n_obs);
  for (i=0; i<n_obs; i++)
    error_mtx[i] = (double *) malloc (sizeof (double) * n_obs);
  ind_seq = (int **) malloc (sizeof (int*) * 2);
  for (i=0; i<2; i++)
    ind_seq[i] = (int *) malloc (sizeof(int)*n_site);
  pi_ind = (double *) malloc (sizeof (double)* n_proc);
  for (i=0; i<n_proc; i++)
    pi_ind[i] = lambda[i];
	
	prob_0to1 = (double **) malloc (sizeof (double *) * n_proc);
	prob_1to0 = (double **) malloc (sizeof (double *) * n_proc);
	for (i=0; i<n_proc; i++)
	{
		prob_0to1[i] = (double *) malloc (sizeof (double) * (n_site-1));
		prob_1to0[i] = (double *) malloc (sizeof (double) * (n_site-1));
	}
		
	for (i=0; i<n_proc; i++)
		for (j=0; j<n_site-1; j++)
		{
			prob_0to1[i][j] = lambda[n_proc + (n_site-1)*i + j];
			prob_1to0[i][j] = lambda[n_proc + (n_site-1)*(n_proc+i) + j];
		}
	
  m = (double *) malloc (sizeof (double) * n_site);
  for (j=0; j<n_site; j++)
    m[j] = lambda[n_proc+2*n_proc*(n_site-1)+j];

	/**************************************
	 * check
	 *************************************/
	if (print_flag)
	{
		printf ("lambda is:\n");
		for (i=0; i<3+3*2*(n_site-1)+n_site; i++)
			printf ("%lf, ", lambda[i]);
		printf ("\n");
		printf ("probs from 0 to 1:\n");
		for (i=0; i<n_proc; i++)
		{
			for (j=0; j<(n_site-1); j++)
				printf ("%lf ", prob_0to1[i][j]);
			printf ("\n");
		}
		printf ("probs from 1 to 0:\n");
		for (i=0; i<n_proc; i++)
		{
			for (j=0; j<(n_site-1); j++)
				printf ("%lf ", prob_1to0[i][j]);
			printf ("\n");
		}
		printf ("m is:\n");
		for (i = 0; i<n_site; i++)
			printf ("%lf ", m[i]);
		printf ("\n");
	}

	/*************************************
	 * compute pi, trans.mtx, emiss.mtx
	 * order: M, RP, RD
	 ************************************/
  /* pi */
  compute_pi (pi, pi_ind);
  if (print_flag)
    {
      printf ("pi_ind:\n");
      for (i=0; i<n_proc; i++)
	printf ("%lf ", pi_ind[i]);
      printf ("\n");
      printf ("pi:\n");
      for (i=0; i<n_hidden; i++)
	printf ("%lf ", pi[i]);
      printf ("\n");
    }

  /* trans.mtx */
  compute_trans_mtx_jump (trans_mtx, prob_0to1, prob_1to0, n_proc, n_site, print_flag);
  if (print_flag)
    {
      printf ("transition matrix:\n");
      for (i=0; i<n_hidden; i++)
	{
	  for (j=0; j<n_hidden; j++)
	    printf ("%lf ", trans_mtx[0][i][j]);
	  printf ("\n");
	}
    }

  /* emission matrix */
  compute_error_mtx (error_mtx, error_b, error_c, print_flag);
  compute_emiss_mtx_c (emiss_mtx, m, m_denovo, m_maint, d_denovo, d_maint, error_mtx, n_site, n_hidden, n_obs, de_novo_flag, print_flag);
  if (print_flag)
    {
		for (s=0; s<n_site; s++)
		{
			printf ("emission matrix (with error) at site %d:\n", s);
			for (i=0; i<n_hidden; i++)
			{
				for (j=0; j<n_obs; j++)
					printf ("%lf ", emiss_mtx[s][i][j]);
				printf ("\n");
			}
		}
    }

  /*************************************
   * compute the log likelihood
   ************************************/
  for (i=0; i<n_seq; i++)
    {
      for (j=0; j<n_site; j++)
	{
	  ind_seq[0][j] = data[i*2*n_site+j];
	  ind_seq[1][j] = data[(i*2+1)*n_site+j];
	}


      if (print_flag)
	{
	  printf ("seq #%d:\n", i);
	  for (k=0; k<2; k++)
	    {
	      for (j=0; j<n_site; j++)
		printf ("%d ", ind_seq[k][j]);
	      printf ("\n");
	    }
	}


      /******************************************
       * compute log likelihood for an unordered 
       * double-stranded sequence
       *****************************************/
      strandtype[i] = loglikind_strandtype_hmm (ind_seq, pi, trans_mtx, emiss_mtx, n_site, n_hidden, print_flag);

    }

	for (s=0; s<n_proc; s++)
	{
		free (prob_0to1[s]);
		free (prob_1to0[s]);
	}
	free (prob_0to1);
	free (prob_1to0);
  for (s=0; s<n_site; s++)
    {
      for (i=0; i<n_hidden; i++)
	  free (emiss_mtx[s][i]);
      free (emiss_mtx[s]);
    }
  for (s=0; s<n_site-1; s++)
    {
      for (i=0; i<n_hidden; i++)
	  free (trans_mtx[s][i]);
      free (trans_mtx[s]);
    }
  free (trans_mtx);
  free (emiss_mtx);
  for (i=0; i<n_obs; i++)
    free (error_mtx[i]);
  free (error_mtx);
  free (pi);
  free (pi_ind);
  for (i=0; i<2; i++)
    free (ind_seq[i]);
  free (ind_seq);
  free (m);
}


/**********************************************************
 * compute_pi
 *
 * function to compute the initial distribution pi
 * from transition probabilities rho's
 *
 * input: pi_ind: initial prob of meth. on each meth. process; 
 *                vector of length n_pi_ind=3
 * output: pi: initial distribution; vector of n_hidden
 **********************************************************/
void compute_pi (double * pi, double * pi_ind)
{
  int i,j,k;

  for (i=0; i<2; i++)
    for (j=0; j<2; j++)
      for (k=0; k<2; k++)
	pi[i*4+j*2+k] = ((1-i)*(1-pi_ind[0])+i*pi_ind[0]) * ((1-j)*(1-pi_ind[1])+j*pi_ind[1]) * ((1-k)*(1-pi_ind[2])+k*pi_ind[2]);

}


/*****************************************************************
 * compute_trans_mtx_jump
 *
 * function to compute transition matrix for the hidden 
 * Markov chain as a jump process; transition probabilties are 
 * different at different sites.
 *
 * input: rho_0: Pr(0->1); matrix of n_proc by n_site-1
 *		  rho_1: Pr(1->0); matrix of n_proc by n_site-1
 *        n_proc: number of methylation processes; equals 3 here
 *		  n_site: number of sites
 *        print_flag: 0 for no output; 1 for output
 * output: trans_mtx: transition matrix; n_site-1 by n_hidden by n_hidden
 *         where n_hidden=8
 *****************************************************************/
void compute_trans_mtx_jump (double *** trans_mtx, double ** rho_0, double ** rho_1, int n_proc, int n_site, int print_flag)
{
  double ** M_trans, ** RP_trans, ** RD_trans, ** R_trans;
  int i,j,t;
  
  M_trans = (double **) malloc (sizeof(double*)*2);
  RP_trans = (double **) malloc (sizeof(double*)*2);
  RD_trans = (double **) malloc (sizeof(double*)*2);
  for (i=0; i<2; i++)
    {
      M_trans[i] = (double *) malloc (sizeof(double)*2);
      RP_trans[i] = (double *) malloc (sizeof(double)*2);
      RD_trans[i] = (double *) malloc (sizeof(double)*2);
    }
	R_trans = (double **) malloc (sizeof(double*)*4);
	for (i=0; i<4; i++)
		R_trans[i] = (double *) malloc (sizeof(double)*4);
	
  for (t=0; t<n_site-1; t++)
  {
	M_trans[0][0] = 1.0-rho_0[0][t];
	M_trans[0][1] = rho_0[0][t];
	M_trans[1][0] = rho_1[0][t];
	M_trans[1][1] = 1.0-rho_1[0][t];
	RP_trans[0][0] = 1.0-rho_0[1][t];
	RP_trans[0][1] = rho_0[1][t];
	RP_trans[1][0] = rho_1[1][t];
	RP_trans[1][1] = 1.0-rho_1[1][t];
	RD_trans[0][0] = 1.0-rho_0[2][t];
	RD_trans[0][1] = rho_0[2][t];
	RD_trans[1][0] = rho_1[2][t];
	RD_trans[1][1] = 1.0-rho_1[2][t];

	Kronecker_prod (RD_trans, 2, RP_trans, 2, R_trans);

	if (print_flag)
	{
		printf ("M_trans:\n");
		for (i=0; i<2; i++)
		{
			for (j=0; j<2; j++)
				printf ("%lf ", M_trans[i][j]);
			printf ("\n");
		}
		printf ("RP_trans:\n");
		for (i=0; i<2; i++)
		{
			for (j=0; j<2; j++)
				printf ("%lf ", RP_trans[i][j]);
			printf ("\n");
		}
		printf ("RD_trans:\n");
		for (i=0; i<2; i++)
		{
			for (j=0; j<2; j++)
				printf ("%lf ", RD_trans[i][j]);
			printf ("\n");
		}
    }

	Kronecker_prod (R_trans, 4, M_trans, 2, trans_mtx[t]);
	  if (print_flag)
	  {
		  printf ("trans:\n");
		  for (i=0; i<8; i++)
		  {
			  for (j=0; j<8; j++)
				  printf ("%lf ", trans_mtx[t][i][j]);
			  printf ("\n");
		  }
	  }		  
  }

  for (i=0; i<2; i++)
  {
      free (M_trans[i]);
      free (RP_trans[i]);
      free (RD_trans[i]);
  }
  free (M_trans);
  free (RP_trans);
  free (RD_trans);
  for (i=0; i<4; i++)
    free (R_trans[i]);
  free (R_trans);
}


void compute_error_mtx (double ** error_mtx, double b, double c, int print_flag)
{
  int i,j;
  double ** error_mtx_tmp;

  error_mtx_tmp = (double **) malloc (sizeof (double *) * 2);
  for (i=0; i<2; i++)
    error_mtx_tmp[i] = (double *) malloc (sizeof (double) * 2);

  error_mtx_tmp[0][0] = 1.0-b;
  error_mtx_tmp[0][1] = b;
  error_mtx_tmp[1][0] = c;
  error_mtx_tmp[1][1] = 1.0-c;
  
  Kronecker_prod (error_mtx_tmp, 2, error_mtx_tmp, 2, error_mtx);

  if (print_flag)
    {
      printf ("2-by-2 error matrix: ");
      for (i=0; i<2; i++)
	for (j=0; j<2; j++)
	  printf ("%lf ", error_mtx_tmp[i][j]);
      printf ("\n");
      printf ("4-by-4 error matrix:\n");
      for (i=0; i<4; i++)
	{
	  for (j=0; j<4; j++)
	    printf ("%lf ", error_mtx[i][j]);
	  printf ("\n");
	}
    }

  for (i=0; i<2; i++)
    free (error_mtx_tmp[i]);
  free (error_mtx_tmp);
}


void compute_emiss_mtx_c (double *** emiss_mtx_c, double * m, double m_denovo, double m_maint, double d_denovo, double d_maint, double ** error_mtx, int n_site, int n_hidden, int n_obs, int DE_NOVO, int print_flag)
{
  double ** emiss_mtx;
  int i,j,s;
  double m_s;
  
  emiss_mtx = (double **) malloc (sizeof (double *) * n_hidden);
  for (i=0; i<n_hidden; i++)
    emiss_mtx[i] = (double *) malloc (sizeof (double) * n_obs);

  /*
  for (s=0; s<n_site; s++)
    for (i=0; i<n_hidden; i++)
      for (j=0; j<n_obs; j++)
	emiss_mtx_c[s][i][j] = 0.0;
  */

  for (i=0; i<n_hidden; i++)
    for (j=0; j<n_obs; j++)
      emiss_mtx[i][j] = 0.0;

  for (s=0; s<n_site; s++)
    {
      m_s = m[s];
      /* assumption 1: de novo events occur only at previously 
	 unmethylated sites */
      if (DE_NOVO==0)
	{
	  emiss_mtx[0][0] = 1.0-m_s;
	  emiss_mtx[0][2] = m_s;
	  emiss_mtx[1][1] = 1.0-m_s;
	  emiss_mtx[1][3] = m_s;
	  emiss_mtx[2][2] = 1.0;
	  emiss_mtx[3][2] = m_s;
	  emiss_mtx[3][3] = 1.0-m_s;
	  emiss_mtx[4][0] = 1.0-m_s;
	  emiss_mtx[4][3] = m_s;
	  emiss_mtx[5][1] = 1.0-m_s;
	  emiss_mtx[5][3] = m_s;
	  emiss_mtx[6][2] = 1.0-m_s;
	  emiss_mtx[6][3] = m_s;
	  emiss_mtx[7][3] = 1.0;
	}

		// assumption 2: de novo activities are allowed for Dnmt1 
		// maintenance and de novo activities of DNMT3 have different rates
      if (DE_NOVO==1)
	{
/*
		emiss_mtx[0][0] = 1.0-m_s;
	  emiss_mtx[0][2] = m_s;
	  emiss_mtx[1][1] = 1.0-m_s;
	  emiss_mtx[1][3] = m_s;
	  emiss_mtx[2][2] = 1.0;
	  emiss_mtx[3][3] = 1.0;
	  emiss_mtx[4][0] = (1.0 - m_s) * (1.0 - m_denovo + m_denovo * m_fail);
	  emiss_mtx[4][1] = (1.0 - m_s) * m_denovo * (1.0 - m_fail);
	  emiss_mtx[4][2] = m_s * m_fail;
	  emiss_mtx[4][3] = m_s * (1.0 - m_fail);
	  emiss_mtx[5][1] = 1.0-m_s;
	  emiss_mtx[5][3] = m_s;
	  emiss_mtx[6][2] = (1.0 - m_s) * (1.0 - m_denovo + m_denovo * m_fail) + m_s * m_fail;
	  emiss_mtx[6][3] = (1.0 - m_s) * m_denovo * (1.0 - m_fail) + m_s * (1.0 - m_fail);
	  emiss_mtx[7][3] = 1.0;
 */
		emiss_mtx[0][0] = 1.0-m_s;
		emiss_mtx[0][2] = m_s;
		emiss_mtx[1][0] = (1.0 - m_s) * (1.0 - d_denovo);
		emiss_mtx[1][1] = (1.0 - m_s) * d_denovo;
		emiss_mtx[1][2] = m_s * (1.0-d_maint);
		emiss_mtx[1][3] = m_s * d_maint;
		emiss_mtx[2][2] = 1.0;
		emiss_mtx[3][2] = (1.0 - m_s) * (1.0 - d_denovo) + m_s * (1.0 - d_maint);
		emiss_mtx[3][3] = (1.0 - m_s) * d_denovo + m_s * d_maint;
		emiss_mtx[4][0] = (1.0 - m_s) * (1.0 - m_denovo);
		emiss_mtx[4][1] = (1.0 - m_s) * m_denovo;
		emiss_mtx[4][2] = m_s * (1.0-m_maint);
		emiss_mtx[4][3] = m_s * m_maint;
		emiss_mtx[5][0] = (1.0-m_s) * (1.0 - d_denovo - d_maint + d_denovo * d_maint);
		emiss_mtx[5][1] = (1.0-m_s) * (d_denovo + d_maint - d_denovo * d_maint);
		emiss_mtx[5][2] = m_s * (1.0 - d_denovo - d_maint + d_denovo * d_maint);
		emiss_mtx[5][3] = m_s * (d_denovo + d_maint - d_denovo * d_maint);
		emiss_mtx[6][2] = (1.0 - m_s) * (1.0 - m_denovo) + m_s * (1.0 - m_maint);
		emiss_mtx[6][3] = (1.0 - m_s) * m_denovo + m_s * m_maint;
		emiss_mtx[7][3] = 1.0;
	}

  if (print_flag)
    {
		printf ("emission matrix (without error):\n");
		for (i=0; i<n_hidden; i++)
		{
			for (j=0; j<n_obs; j++)
				printf ("%lf ", emiss_mtx[i][j]);
			printf ("\n");
		}
	}

      /* multiply emiss_mtx and the error matrix */
      multiply_mtx (emiss_mtx_c[s], emiss_mtx, error_mtx, n_hidden, n_obs, n_obs, print_flag);
    }

  for (i=0; i<n_hidden; i++)
    free (emiss_mtx[i]);
  free (emiss_mtx);
}


/***********************************************************************
 * multiply_mtx
 * 
 * function compute the product of two matrices of different dimensions
 *
 * input: mtx_1: left matrix; dim_i by dim_j
 *        mtx_2: right matrix; dim_j by dim_k
 *        dim_i: number of rows of mtx_1
 *        dim_j: number of columns of mtx_1, or number of rows of mtx_2
 *        dim_k: number of columns of mtx_2
 *        int print_flag: 1 for output and 0 for no output
 * output: mtx_prod: product of mtx_1 and mtx_2; dim_i by dim_k
 ***********************************************************************/
void multiply_mtx (double ** mtx_prod, double ** mtx_1, double ** mtx_2, int dim_i, int dim_j, int dim_k, int print_flag)
{
  int i,j,k;

  for (i=0; i<dim_i; i++)
    for (k=0; k<dim_k; k++)
      {
		mtx_prod[i][k] = 0.0;
		for (j=0; j<dim_j; j++)
			mtx_prod[i][k] += mtx_1[i][j] * mtx_2[j][k];
      }

  if (print_flag)
    {
      printf ("mtx_prod:\n");
      for (i=0; i<dim_i; i++)
		{
			for (k=0; k<dim_k; k++)
				printf ("%lf ", mtx_prod[i][k]);
			printf ("\n");
		}
    }
}


/*********************************************
 * Kronecker_prod
 *
 * Function to compute the Kronecker product
 * of two matrices
 ********************************************/
void Kronecker_prod (double ** mtx1, int dim1, double ** mtx2, int dim2, double ** mtx)
{
  int i, j, m, n;

  for (i = 0; i < dim2; i++)
    for (j = 0; j < dim2; j++)
      for (m = 0; m < dim1; m++)
        for (n = 0; n < dim1; n++)
          mtx[i * dim1 + m][j * dim1 + n] = mtx1[m][n] * mtx2[i][j];
}


/********************************************
 * loglikind_strandtype_hmm
 *
 * function to compute strand type prob of an 
 * unordered double-stranded sequence with 
 * transition (A) and emission (B) probabilities 
 * different at different sites.
 *
 * input: ind_seq: an unordered ds seq; matrix of 2 by length
 *        pi: initial distribution; vector of length "length"
 *        A: transition matrix; length-1 by n_hidden by n_hidden
 *        B: emission matrix; length by n_hidden by nObsState
 *        length: length of observed sequence
 *        n_hidden: number of hidden states
 *        print_flag
 * output: ind_strand_type: prob of top strand being parent
 *******************************************/
double loglikind_strandtype_hmm (int ** ind_seq, double * pi, double *** A, double *** B, int length, int n_hidden, int print_flag)
{
  int flag_tb = 0;  /* flag of whether top and bottom strands differ */
  int ** alt_seq;
  double likind;
  int i,j;
  double prob_top_parent, ind_strand_type;

  if (print_flag)
    {
      for (i=0; i<2; i++)
	{
	  for (j=0; j<length; j++)
	    printf ("%d ", ind_seq[i][j]);
	  printf ("\n");
	}
    }

  alt_seq = (int **) malloc (sizeof (int *) * 2);
  for (i=0; i<2; i++)
    alt_seq[i] = (int *) malloc (sizeof (int) * length);

  /**********************************************
   * compute likelihood of the sequence with 
   * top parent and bottom daughter
   *********************************************/
  likind = likind_ordered_hmm (ind_seq, pi, A, B, length, n_hidden, print_flag);
  prob_top_parent = likind;

  /***************************************************
   * compare top and bottom strands
   * if they differ, compute likelihood of the 
   * other orientation (top daughter, bottom parent)
   **************************************************/
  for (j=0; j<length; j++)
    {
      if (ind_seq[0][j] != ind_seq[1][j])
		flag_tb += 1;
    }
  if (flag_tb>0)          /* flag_tb>0 if top!=bottom */
    {
      /*************************************
       * swap top and bottom strands
       ************************************/
      for (i=0; i<2; i++)
		for (j=0; j<length; j++)
			alt_seq[i][j] = ind_seq[1-i][j];
      likind += likind_ordered_hmm (alt_seq, pi, A, B, length, n_hidden, print_flag);
    }
	
  if (flag_tb>0)
  {
	ind_strand_type = prob_top_parent / likind;
	  if (print_flag)
	  {
		  printf ("prob_top_parent: %lf\n", prob_top_parent);
		  printf ("lik ind: %lf\n", likind);
		  printf ("ind_strand_type: %lf\n", ind_strand_type);
	  }
  }
  else
    ind_strand_type = 0.5;

  for (i=0; i<2; i++)
    free (alt_seq[i]);
  free (alt_seq);

  return (ind_strand_type);
}


/********************************************************************
 * likind_ordered_hmm
 *
 * function to compute likelihood of an ordered double-stranded
 * sequence with transition (A) and emission (B) probabilities different 
 * at different sites.
 *
 * input: ind_seq: an ordered double-stranded sequence; matrix of 2 by length
 *        pi: initial distribution; vector of length nHiddenState
 *        A: transition matrix; matrix of length-1 by nHiddenState by nHiddenState
 *        B: emission matrix; matrix of length by nHiddenState by nObsState
 *        length: length of the observed sequence
 *        nHiddenState: number of hidden states
 *        print_flag
 * output: likind: scalar; likelihood of an ordered double-stranded 
 *                 sequence
 *******************************************************************/
double likind_ordered_hmm (int ** ind_seq, double * pi, double *** A, double *** B, int length, int nHiddenState, int print_flag)
{
  int j;
  double likind;
  int * Obs;

  Obs = (int *) malloc (sizeof(int) * length);

  /**************************************
   * re-code the double-stranded sequence
   * (0,0)-0, (0,1)-1, (1,0)-2, (1,1)-3
   *************************************/
  for (j=0; j<length; j++)
    Obs[j] = ind_seq[0][j]*2+ind_seq[1][j];

  /****************************************
   * compute the likelihood
   ***************************************/
  likind = forward_backward (pi, A, B, Obs, length, nHiddenState, print_flag);

  free (Obs);
  return (likind);

}


/*********************************************************************
 * forward_backward
 * 
 * function to compute likelihood for an observed sequence under a 
 * hidden Markov model with transition and emission probabilities 
 * different at different sites.
 *
 * input: pi: initial distribution; vector of length nHiddenState
 *        A: transition matrix; matrix of length-1 by nHiddenState by nHiddenState
 *        B: emission matrix; matrix of length by nHiddenState by nObsState
 *        Obs: observed sequence; vector of length "length"
 *        length: length of observed sequence
 *        nHiddenState: number of hidden states
 *        print_flag: 1 for output and 0 o.w.
 * output: L: likelihood of the sequence
 *********************************************************************/
double forward_backward (double * pi, double *** A, double *** B, int * Obs, int length, int nHiddenState, int print_flag)
{ 

  int t, i, j, k;
  double ** Alpha;            /* forward pass */
  double ** Beta;             /* backward pass */
  double L;

  if (print_flag)
  {
    printf ("length:%d\n", length);
	printf ("observed sequence:\n");
	for (i=0; i<length; i++)
		printf ("%d ", Obs[i]);
	printf ("\n");
  }

  Alpha = (double **) malloc (sizeof (double *) * length);
  for (t = 0; t < length; t++)
    Alpha[t] = (double *) malloc (sizeof (double) * nHiddenState);
  Beta = (double **) malloc (sizeof (double *) * length);
  for (t = 0; t < length; t++)
    Beta[t] = (double *) malloc (sizeof (double) * nHiddenState);

  /*******************************************
   * calcualte the forward pass 
   ******************************************/
  for (t = 0; t < length; t++)
	{
		if (t == 0)
		{
			for (i = 0; i < nHiddenState; i++)
			{
				Alpha[t][i] = pi[i] * B[t][i][Obs[t]];
			}
		}
		else
		{
			for (i = 0; i < nHiddenState; i++)
			{
				Alpha[t][i] = 0.0;
				for (j = 0; j < nHiddenState; j++)
				{
					Alpha[t][i] += Alpha[t - 1][j] * A[t-1][j][i] * B[t][i][Obs[t]];
				}
			}
		}
    }

  if (print_flag)
    {
		printf ("log of alpha:\n");
		for (t = 0; t < length; t++)
		{
			printf ("site %d: ", t);
			for (i = 0; i < nHiddenState; i++)
				printf ("%4.3lf ", log (Alpha[t][i]));
			printf ("\n");
		}
    }

  /*****************************************
   * calculate the backward pass 
   ****************************************/
  for (t = length - 1; t >= 0; t--)
    {
		if (t == length - 1)
			for (i = 0; i < nHiddenState; i++)
				Beta[t][i] = 1.0;
		else
		{
			for (i = 0; i < nHiddenState; i++)
			{
				Beta[t][i] = 0.0;
				for (j = 0; j < nHiddenState; j++)
				{
					Beta[t][i] += A[t][i][j] * B[t+1][j][Obs[t + 1]] * Beta[t + 1][j];
				}
			}
		}
    }


  /*******************************************
   * calculate the total likelihood 
   ******************************************/
  L = 0.0;
  for (i = 0; i < nHiddenState; i++)
    L += Alpha[length - 1][i];
  if (print_flag)
    printf ("log likelihood is: %lf\n", log(L));

  if (print_flag)
    {
		printf ("log likelihood at each site:\n");
		for (t = 0; t < length; t++)
		{
			L = 0.0;
			for (i = 0; i < nHiddenState; i++)
				L += Alpha[t][i] * Beta[t][i];
			printf ("site %d: %lf\n", t, log(L));
		}
    }

  for (t = 0; t < length; t++)
    {
      free (Alpha[t]);
      free (Beta[t]);
    }
  free (Alpha);
  free (Beta);

  return (L);
}

