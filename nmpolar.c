/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
NMPOLAR.C - Functions for polarization analysis of three-component data 
*************************************************************************

Author: Nils Maercklin,
    RISSC, University of Naples, Italy, September 2009

Version: 2010-03-16

Modifications:
    2009-09-15 (NM): first version based on functions from SUPOLAR
    2010-01-04 (NM): added normalized inclination angles calc_norminc()
    2010-03-16 (NM): eigenanalysis function do_eigen_zm() for zero-mean data

Most of these functions are taken from SUPOLAR (with modifications). 
Polarization attribute definitions are given in Maercklin (1999) and 
in a SUPOLAR manual available at http://purl.org/net/nils/man/supolar.

Note that all array indices are zero-based here, whereas corresponding 
indices in SUPOLAR may be one-based.

The functions do_eigen() and do_eigen_zm() require memory allocation 
functions defined in "nmutils.c". All other functions use only standard 
C libraries.


References:
Benhama, A., Cliet, C., and Dubesset, M. (1988). Study and applications of 
    spatial directional filterings in three-component recordings. 
    Geophys. Prosp., 36(6), 591-613.
Jurkevics, A. (1988). Polarization analysis of three-component array data. 
    Bull. Seism. Soc. Am., 78(5), 1725-1743.
Kanasewich, E. R. (1981). Time Sequence Analysis in Geophysics.
    The University of Alberta Press.
Maercklin, N. (1999). Polarisationsanalyse refraktionsseismischer Daten vom 
    Vulkan Merapi, Indonesien. Masters thesis (Diplomarbeit), Geophysics, 
    Christian-Albrechts-University, Kiel, Germany.
Montalbetti, J. R. and Kanasewich, E. R. (1970). Enhancement of 
    teleseismic body phases with a polarization filter. 
    Geophys. J. R. Astr. Soc., 21(2), 119-129.
Press, W. H., Teukolsky, S. A., Vetterling, W. T., and Flannery, B. P. 
    (1996). Numerical Recipes in C - The Art of Scientific Computing.
    Cambridge University Press, Cambridge.
Samson, J. C. (1973). Descriptions of the polarization states of vector 
    processes: applications to ULF magnetic fields. 
    Geophys. J. R. Astr. Soc., 34(4), 403-419.
*************************************************************************/

/* Include files */
#include <stdlib.h>
#include <math.h>
#include "nmutils.h"


/* Function prototypes */
void eig_jacobi(float **a, float d[], float **v, int n);
void sort_eigenvalues(float d[], float **v, int n);

float covar(float *data1, float *data2, int istart, int iwl);
void do_eigen(float **data, int nt, int iwl, float **ddata, float **vdata);
void do_eigen_zm(float **data, int nt, int iwl, float **ddata, float **vdata);

float calc_rl(float *d, float rlq, int opt);
float calc_ellip(float *d, int i1, int i2);
float calc_f1(float *d);
float calc_l1(float *d);
float calc_plan(float *d);
float calc_tau(float *d);
float calc_er(float *d);
float calc_norminc(float *v);
float calc_theta(float *v, int opt);
float calc_phi(float *v, int opt);

#ifndef PI
#define PI (3.141592653589793)
#endif


/************************************************************************
Functions for polarization analysis based on the covariance matrix
*************************************************************************/

float covar(float *data1, float *data2, int istart, int iwl)
/************************************************************************
covar - covariance in a boxcar time window

Input:
data1   time series, array of floats
data2   time series, array of floats
istart  index of first sample in analysis time window
iwl     analysis time window length in samples

Output: returns covariance

Reference:
Montalbetti, J. R. and Kanasewich, E. R. (1970). Enhancement of 
    teleseismic body phases with a polarization filter. 
    Geophys. J. R. Astr. Soc., 21(2), 119-129.

Author: Nils Maercklin, 1998
*************************************************************************/
{
    register int i;
    float cov=0.0;
    float mean1=0.0;
    float mean2=0.0;
    
    for (i=istart; i<(istart+iwl); i++) {
        mean1 += data1[i];
        mean2 += data2[i];
    }
    mean1 = mean1 / (float) iwl;
    mean2 = mean2 / (float) iwl;

    for (i=istart; i<(istart+iwl); i++) {
        cov += (data1[i]-mean1) * (data2[i]-mean2);
    }
    cov = cov / ((float) iwl);

    return cov;
}



void do_eigen(float **data, int nt, int iwl, float **ddata, float **vdata)
/************************************************************************
do_eigen - compute eigenvalues and eigenvectors in moving time window

Input:
data    three-component data, array of floats (data[3][nt])
nt      number of samples per trace
iwl     analysis time window length in samples

Output:
ddata   eigenvalues, array of three-component float vectors (ddata[nt][3])
vdata   eigenvectors of largest eigenvalue ddata[i][0]  (vdata[nt][9])

Author: Nils Maercklin, 1998 - 2009
*************************************************************************/
{
    register int i, j, it;
    float **a=NULL;    /* covariance matrix */
    float **v=NULL;    /* matrix of eigenvectors */
    float  *d=NULL;    /* array of eigenvalues */


    /* Allocate space */
    a = fmalloc2(3,3);
    v = fmalloc2(3,3);
    d = fmalloc1(3);

    /* Loop over samples */
    for (it=iwl;it<nt;it++) {
        /* Covariance matrix */
        for (i=0;i<3;i++) {
            for (j=i;j<3;j++) {
                a[i][j] = a[j][i] = \
                    covar(data[i], data[j], it-iwl, iwl);
            }
        }

        /* Compute eigenvalues and -vectors */
        eig_jacobi(a, d, v, 3);
        sort_eigenvalues(d, v, 3);

        /* Save eigenvalues and eigenvectors */
        for (i=0;i<3;i++) {
            ddata[it][i]   = d[i];
            vdata[it][i]   = v[i][0];
            vdata[it][i+3] = v[i][1];
            vdata[it][i+6] = v[i][2];
        }
    }

    /* Fill */
    for (it=0;it<iwl;it++) {
        for (i=0;i<3;i++) ddata[it][i] = ddata[iwl][i];
        for (i=0;i<9;i++) vdata[it][i] = vdata[iwl][i];
    }

    /* Free space */
    free2((void *)a);
    free2((void *)v);
    free(d);
}



void do_eigen_zm(float **data, int nt, int iwl, float **ddata, float **vdata)
/************************************************************************
do_eigen_zm - compute eigenvalues and eigenvectors in moving time window
              efficiently for zero-mean data

Input:
data    three-component data, array of floats (data[3][nt]), zero mean
nt      number of samples per trace
iwl     analysis time window length in samples

Output:
ddata   eigenvalues, array of three-component float vectors (ddata[nt][3])
vdata   eigenvectors of largest eigenvalue ddata[i][0]  (vdata[nt][9])

Note:
This function assumes that the data in all analysis windows of length iwl 
are zero-mean. This function runs much faster than do_eigen().

Author: Nils Maercklin, 2010
*************************************************************************/
{
    register int i, j, it;
    float **a=NULL;    /* covariance matrix */
    float **v=NULL;    /* matrix of eigenvectors */
    float  *d=NULL;    /* array of eigenvalues */
    float **s=NULL;    /* accumulator for covariance matrix */
    float c=0.0;       /* scale factor */

    /* Allocate space */
    a = fmalloc2(3,3);
    v = fmalloc2(3,3);
    d = fmalloc1(3);
    s = fmalloc2(3,3);

    /* Initialize */
    c = (iwl) ? 1.0 / ((float)iwl) : 1.0;
    for (i=0; i<3; i++) {
        for (j=i; j<3; j++) s[i][j] = s[j][i] = 0.0;
    }

    /* Loop over samples */
    for (it=0; it<nt; it++) {
        /* Update covariance matrix */
        if (it >= iwl) {
            for (i=0; i<3; i++) {
                for (j=i; j<3; j++) {
                    s[i][j] -= c * data[i][it-iwl] * data[j][it-iwl];
                }
            }
        }
        for (i=0; i<3; i++) {
            for (j=i; j<3; j++) {
                s[i][j] += c * data[i][it] * data[j][it];
                a[i][j] = a[j][i] = s[i][j];
            }
        }

        /* Compute eigenvalues and -vectors */
        eig_jacobi(a, d, v, 3);
        sort_eigenvalues(d, v, 3);

        /* Save eigenvalues and eigenvectors */
        for (i=0;i<3;i++) {
            ddata[it][i]   = d[i];
            vdata[it][i]   = v[i][0];
            vdata[it][i+3] = v[i][1];
            vdata[it][i+6] = v[i][2];
        }
    }

    /* Free space */
    free2((void *)a);
    free2((void *)v);
    free2((void *)s);
    free(d);
}



float calc_rl(float *d, float rlq, int opt)
/************************************************************************
calc_rl - compute rectilinearity RL (polarization attribute)

Input:
d       three-element array of eigenvalues (d[3])
rlq     contrast parameter 
opt     flag for different definitions of RL

References:
Jurkevics, A. (1988). Polarization analysis of three-component array data. 
    Bull. Seism. Soc. Am., 78(5), 1725-1743.
Kanasewich, E. R. (1981). Time Sequence Analysis in Geophysics.
    The University of Alberta Press.

Author: Nils Maercklin, 1998 - 2009
*************************************************************************/
{
    float rl;

    if (d[0]) {
        if (opt) {
            /* RL definition after Jurkevics (1988) */
            rl = 1.0 - pow( 0.5*(fabs(d[1]/d[0]) + fabs(d[2]/d[0]) ), rlq);
        }
        else {
            /* RL definition after Kanasewich (1981) */
            rl = 1.0 - pow(fabs(d[1]/d[0]), rlq);
        }
        return rl;
    }
    else {
        return 0.0;
    }
}



float calc_ellip(float *d, int i1, int i2)
/************************************************************************
calc_ellip - compute ellipticities (polarization attributes)

Input:
d       three-element array of eigenvalues (d[3])
i1      eigenvalue index, numerator (values 0,1,2)
i2      eigenvalue index, denominator (values 0,1,2)

Author: Nils Maercklin, 1998
*************************************************************************/
{
    float ellip;

    if (d[i2]) {
        ellip = sqrt( fabs(d[i1]/d[i2]) );
        return ellip;
    }
    else
        return 0.0;
}



float calc_f1(float *d)
/************************************************************************
calc_f1 - compute flatness coefficient f1 (polarization attribute)

Input:
d       three-element array of eigenvalues (d[3])

Reference:
Benhama, A., Cliet, C., and Dubesset, M. (1988). Study and applications 
    of spatial directional filterings in three-component recordings. 
    Geophys. Prosp., 36(6), 591-613.

Author: Nils Maercklin, 1998
*************************************************************************/
{
    float f1,x1,x2;

    x1 = 3.0 * calc_ellip(d,2,0);
    x2 = 1.0 + calc_ellip(d,1,0) + calc_ellip(d,2,0);
    f1 = 1.0 - x1 / x2; 
    return f1;
}




float calc_l1(float *d)
/************************************************************************
calc_l1 - compute linearity coefficient l1 (polarization attribute)

Input:
d       three-element array of eigenvalues (d[3])

Author: Nils Maercklin, 1998
*************************************************************************/
{
    float l1,x1,x2;

    x1 = 3.0 * ( calc_ellip(d,1,0) + calc_ellip(d,2,0) );
    x2 = 2.0 * ( 1.0 + calc_ellip(d,1,0) + calc_ellip(d,2,0));
    l1 = 1.0 - x1 / x2;
    return l1;
}



float calc_plan(float *d)
/************************************************************************
calc_plan - compute planarity (polarization attribute)

Input:
d       three-element array of eigenvalues (d[3])

Author: Nils Maercklin, 2001
*************************************************************************/
{
    float pln;

    if (d[0]+d[1]) {
        pln = 1.0 - 2.0*d[2] / (d[0] + d[1]);
        return pln;
    }
    else
        return 0.0;
}



float calc_tau(float *d)
/************************************************************************
calc_tau - compute global polarization parameter tau (polarization attribute)

Input:
d       three-element array of eigenvalues (d[3])

Reference:
Samson, J. C. (1973). Descriptions of the polarization states of vector 
    processes: applications to ULF magnetic fields. 
    Geophys. J. R. Astr. Soc., 34(4), 403-419.

Author: Nils Maercklin, 1998
*************************************************************************/
{
    float x1, x2, x3, x4, tau;

    if (d[0]) {
        x1  = (d[0] - d[1]) * (d[0] - d[1]);
        x2  = (d[0] - d[2]) * (d[0] - d[2]);
        x3  = (d[1] - d[2]) * (d[1] - d[2]);
        x4  = (d[0] + d[1] + d[2]) * (d[0] + d[1] + d[2]);
        tau = sqrt( (x1 + x2 + x3) / (2.0*x4) );

        return tau;
    }
    else {
        return 0.0;
    }
}



float calc_er(float *d)
/************************************************************************
calc_er - compute eigenresultant (amplitude in principal polarization dir.)

Input:
d       three-element array of eigenvalues (d[3])

Author: Nils Maercklin, 1998
*************************************************************************/
{
    float er;

    er = sqrt(fabs(d[0]));
    return er;
}



float calc_norminc(float *v)
/************************************************************************
calc_norminc - calculate normalized inclination angles

Input:
v       three-element vector Vi (v[3]), e.g. an eigenvector

Reference:
Jurkevics, A. (1988). Polarization analysis of three-component array data. 
    Bull. Seism. Soc. Am., 78(5), 1725-1743.

Author: Nils Maercklin, 2009
*************************************************************************/
{
    float l, val;

    l = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

    if (l) {
        val = 2.0 * acos( fabs(v[0])) / PI;
    }
    else {
        val = 0.0;
    }

    return val;
}



float calc_theta(float *v, int opt)
/************************************************************************
calc_theta - compute polarization angle theta (incidence angle) in radians

Input:
v       three-element array, eigenvector of principal polarization direction
opt     flag for different definitions of RL

Note:
The vertical component must be v[0]. Output ranges of values are 
[0, pi/2] by default and [-pi/2, pi/2] for opt=1.

References:
Jurkevics, A. (1988). Polarization analysis of three-component array data. 
    Bull. Seism. Soc. Am., 78(5), 1725-1743.
Kanasewich, E. R. (1981). Time Sequence Analysis in Geophysics.
    The University of Alberta Press.

Author: Nils Maercklin, 1998 - 2009
*************************************************************************/
{
    float theta, horiz;

    if (opt==1) {
        /* Definition after Kanasewich (1981) */
        if (v[0]) {
            horiz = sqrt( v[1]*v[1] + v[2]*v[2] );
            theta = atan( horiz / v[0] );
        }
        else {
            theta = 0.0;
        }
    }
    else {
        /* Definition after Jurkevics (1988) */
        theta = acos( fabs(v[0]) );
    }

    return theta;
}



#define VSIGN ( (v[0]<0.0) ? -1.0 : 1.0 )
float calc_phi(float *v, int opt)
/************************************************************************
calc_theta - compute polarization angle phi (azimuth) in radians

Input:
v       three-element array, eigenvector of principal polarization direction
opt     flag for different definitions of RL

Note:
The vertical component must be v[0]. Output ranges of values are 
[-pi/2, pi/2] for opt=1, [0, 2*pi] for opt=3, and [-pi, pi] by default.
The SIGN function is introduced to resolve the 180 deg ambiguity by 
taking the positive vertical component of v[0] (Jurkevics, 1988).

References:
Jurkevics, A. (1988). Polarization analysis of three-component array data. 
    Bull. Seism. Soc. Am., 78(5), 1725-1743.
Kanasewich, E. R. (1981). Time Sequence Analysis in Geophysics.
    The University of Alberta Press.

Author: Nils Maercklin, 1998 - 2009
*************************************************************************/
{
    float phi;

    if (opt==1) {
        /* Definition after Kanasewich (1981), -pi/2 ... pi/2 */
        if (v[1]) {
            phi = atan( v[2] / v[1] );
        }
        else {
            phi = (v[2]>0.0) ? 0.5*PI : -0.5*PI;
        }
    }
    else {
        /* Definition after Jurkevics (1988), -pi ... pi */
        if (v[1]) {
            phi = atan2( v[2]*VSIGN, v[1]*VSIGN);
        }
        else {
            phi = (v[2]>0.0) ? 0.5*PI*VSIGN : -0.5*PI*VSIGN;
        }

        /* Definition after Jurkevics (1988), 0.0 ... 2*pi */
        if (phi<0.0 && opt==3) phi += 2.0*PI;
    }

    return phi;
}
#undef VSIGN



/************************************************************************/
/* Eigenvalues and eigenvectors                                         */
/************************************************************************/

/* Macro used internally (function "eig_jacobi()" */
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
    a[k][l]=h+s*(g-h*tau);

void eig_jacobi(float **a, float d[], float **v, int n)
/**********************************************************************
eig_jacobi - find eigenvalues and corresponding eigenvectors via 
             the jacobi algorithm for symmetric matrices

Input:
a       symmetric matrix a[n][n], e.g. covariance matrix a[3][3]
n       dimension of matrix

Output:
d       array d[n] of eigenvalues
v       array v[n][n] of corresponding eigenvectors

Note:   Array indices are in the range [0, n-1].

Reference:
Press, W. H., Teukolsky, S. A., Vetterling, W. T., and Flannery, B. P. 
    (1996). Numerical Recipes in C - The Art of Scientific Computing.
    Cambridge University Press, Cambridge.

Credits: inspired by Press et al., modifications for CWP/SU and by NM
**********************************************************************/
{
    int j,iq,ip,i;
    float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

    /* allocate space temporarily */
    b = malloc(n*sizeof(float));
    z = malloc(n*sizeof(float));

    /* initialize v to the identity matrix */
    for (ip=0;ip<n;ip++) {
        for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
        v[ip][ip]=1.0;
    }
    /* initialize to the diagonal on matrix a */
    for (ip=0;ip<n;ip++) {
        b[ip]=d[ip]=a[ip][ip];
        z[ip]=0.0;
    }

    /* main iteration loop */
    for (i=1;i<=50;i++) {

        sm=0.0;
        for (ip=0;ip<n-1;ip++) {
            for (iq=ip+1;iq<n;iq++)
                sm += fabs(a[ip][iq]);
        }
        /* normal return */
        if (sm == 0.0) {
            free(z);
            free(b);
            return;
        }
        /* tresh values for first 3 sweeps and therafter */
        if (i < 4)
            tresh=0.2*sm/(n*n);
        else
            tresh=0.0;
        for (ip=0;ip<n-1;ip++) {
            for (iq=ip+1;iq<n;iq++) {
                g=100.0*fabs(a[ip][iq]);
                if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
                    && (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
                    a[ip][iq]=0.0;
                else if (fabs(a[ip][iq]) > tresh) {
                    h=d[iq]-d[ip];
                    if ((float)(fabs(h)+g) == (float)fabs(h))
                        t=(a[ip][iq])/h;
                    else {
                        theta=0.5*h/(a[ip][iq]);
                        t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                        if (theta < 0.0) t = -t;
                    }
                    c=1.0/sqrt(1+t*t);
                    s=t*c;
                    tau=s/(1.0+c);
                    h=t*a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq]=0.0;

                    /* Jacobi rotations */
                    for (j=0;j<ip-1;j++) {
                        ROTATE(a,j,ip,j,iq)
                    }
                    for (j=ip+1;j<=iq-1;j++) {
                        ROTATE(a,ip,j,j,iq)
                    }
                    for (j=iq+1;j<n;j++) {
                        ROTATE(a,ip,j,iq,j)
                    }
                    for (j=0;j<n;j++) {
                        ROTATE(v,j,ip,j,iq)
                    }
                }
            }
        }
        for (ip=0;ip<n;ip++) {
            b[ip] += z[ip];
            d[ip]=b[ip];
            z[ip]=0.0;
        }
    }

    /* this will not happen, hopefully */
    fprintf(stderr,"%s: jacobi iteration does not converge\n", __FILE__);
}
#undef ROTATE



void sort_eigenvalues(float d[], float **v, int n)
/**********************************************************************
sort_eigenvalues - sort eigenvalues and corresponding eigenvectors
                   in descending order

Input:
d       array d[n] of eigenvalues (from "eig_jacobi()")
v       array v[n][n] of corresponding eigenvectors (from "eig_jacobi()")

Output: modified d and v, sorted in descending order

Note:   Array indices are in the range [0, n-1].

Reference:
Press, W. H., Teukolsky, S. A., Vetterling, W. T., and Flannery, B. P. 
    (1996). Numerical Recipes in C - The Art of Scientific Computing.
    Cambridge University Press, Cambridge.

Credits: inspired by Press et al., modifications for CWP/SU and by NM
**********************************************************************/
{
    int k,j,i;
    float p;

    for (i=0;i<n-1;i++) {
        p=d[k=i];
        for (j=i+1;j<n;j++)
            if (d[j] >= p) p=d[k=j];
        if (k != i) {
            d[k]=d[i];
            d[i]=p;
            for (j=0;j<n;j++) {
                p=v[j][i];
                v[j][i]=v[j][k];
                v[j][k]=p;
            }
        }
    }
}

/* END OF FILE */
