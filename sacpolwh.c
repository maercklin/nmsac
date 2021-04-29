/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
SACPOLWH - SAC Test program for polarization attributes defined by 
           Wu & Horiuchi (2008) for S-wave identification
           (based on SACPOLAR v.2009-09-11)
*************************************************************************

Author: Nils Maercklin,
    RISSC, University of Naples, Italy, February 2010

Version: 2011-06-12  (experimental code)

Modifications:
  2010-02-08 (NM): First experimental release of this code
  2010-03-17 (NM): Added option -z for fast eigenanalysis of zero-mean 
                   data, and enabled optional zerophase bandpass
  2011-06-12 (NM): Added user-specified output directory

Notes: Input files must be evenly-sampled time series files in SAC 
    binary format (assuming NVHDR=6, IFTYPE=ITIME (1), LEVEN=1).
    Three subsequent files/traces are considered as one 3-C dataset with 
    common trace length, start time, and sampling rate.

    No station/component consistency checks are made (yet), and the code 
    may contain some (undocumented) experimental features.

Reference:
Wu, C. and Horiuchi, S. (2008). Automatic determination of source 
    parameters of the 2007 Noto Hanto earthquake. Earth, Planets, and Space, 
    60, 1053-1057.
*************************************************************************/

#include "nmsaclib.h"     /* SAC-related functions, and other stuff */


/* Self-documentation */
const char *sdoc[] = {
"                                                                       ",
" SACPOLWH - SAC Polarization analysis of 3-C data (Wu & Horiuchi, 2008)",
"                                                                       ",
" Usage: sacpolwh [-p] attribute [...] [parameters] -f sac_files        ",
"        sacpolwh <stdin [-p] attribute [...] [parameters] > stdout     ",
"                                                                       ",
" Parameters and defaults:                                              ",
"   -f        list of SAC binary waveform files (or read from stdin)    ",
"   -p        list of one or more of the polarization attributes        ",
"             kp, k1, k2, k3, and ki (ki = k1^2 * k2^2 * k3^2)          ",
"   -w   0.5  correlation time window length in seconds                 ",
"   -z        flag: assume zero mean in correlation windows (faster)    ",
"   -v        flag: verbose operation                                   ",
"   -dir      output directory (default is same as input directory)     ",
"                                                                       ",
"   -b1    0  Butterworth low-cut frequency in Hz  (0 = no filter)      ",
"   -b2    0  Butterworth high-cut frequency in Hz (0 = no filter)      ",
"   -bp    3  number of poles of Butterworth filters                    ",
"   -bz       flag: apply zerophase (two-pass) Butterworth filter       ",
"                                                                       ",
" The polarization attributes are computed in a moving time window from ",
" the eigenvalues and the principal eigenvector of the three-component  ",
" covariance matrix.                                                    ",
" Three subsequent files/traces are considered as one three-component   ",
" dataset with common trace length, start time, and sampling rate.      ",
" The vertical component must be the first trace, and it must contain   ",
" a P-wave arrival pick in the SAC header field A (attributes k1, k3).  ",
" Each output file has a suffix corresponding to the attribute name.    ",
"                                                                       ",
" NM, 2011-06-12  (beta)                                                ",
"                                                                       ",
NULL};


/* Prototype of a functions used internally */
float dot_product(float *v1, float *v2);
float calc_wuho_kp(float *d, float a);
float calc_wuho_k1(float *v, float *vp);
float calc_wuho_k2(float *d, float a);
float calc_wuho_k3(float **data3, float *vp, int it0, int iwl);
int search_pvector1(float **vdata, float **ddata, int nt, int iwl, int itp, \
    float a, float *vp);

/* Definitions */
#define RAD2DEG  57.29577951

int 
main(int argc, char **argv)
{
    /* Variables */
    register int i,j,it,iarg;/* loop indices */
    SACHEAD hd;          /* SAC header */
    SACHEAD hd3[3];      /* SAC header, three components */
    float *data=NULL;    /* SAC input trace data */
    float *odata=NULL;   /* SAC output trace data */
    float **data3=NULL;  /* SAC 3-C input data array */
    float **vdata=NULL;  /* eigenvector array */
    float **ddata=NULL;  /* eigenvalue array */
    char *outfname=NULL; /* output file name */
    char *outdir=NULL;   /* optional output directory name */
    int isfile=0;        /* file names flag */
    int isvalid=0;       /* internal validation flag (pol. attribute) */
    int swap;            /* byte-swapping flag */
    int verbose=0;       /* verbose flag */
    int zeromean=0;      /* zero-mean flag (assume data have zero mean) */
    int nstat=0;         /* 3-C station counter (for user info) */
    int icomp=0;         /* component counter */

    int nt=0;            /* number of time samples */
    float dt;            /* time sampling interval in seconds */
    float fmin, fmax;    /* Butterworth filter frequencies */
    int npoles;          /* number of poles of Butterworth filters */
    int zerophase;       /* flag for zerophase Butterworth filters */

    float wl=0.0;        /* averaging time window length in seconds */
    int iwl=0;           /* time window length in samples */

    int pstart=1;        /* start of polarization attribute list */
    float rlq=0.0;       /* contrast factor of rectilinearity RL */

    int itp=0;           /* sample index of P-pick */
    float a;             /* adjusting parameter a */
    float vp[3];         /* P-wave polarization vector */


    /* Print documentation */
    if (argc==1) {
        i=0;
        while (sdoc[i]) fprintf(stderr, "%s\n", sdoc[i++]);
        return (EXIT_FAILURE);
    }

    /* Command line parameters */
    if (!fgetpar(argc, argv, "-w",  &wl))     wl=0.5;
    if (!fgetpar(argc, argv, "-q",  &rlq))    rlq=1.0;

    if (!fgetpar(argc, argv, "-b1", &fmin))   fmin=0.0;
    if (!fgetpar(argc, argv, "-b2", &fmax))   fmax=0.0;
    if (!igetpar(argc, argv, "-bp", &npoles)) npoles=3;
    if (!(zerophase=getflag(argc, argv, "-bz"))) zerophase=0;

    if (!(zeromean=getflag(argc, argv, "-z")) && \
        !(zeromean=getflag(argc, argv, "-zm")))  zeromean=0;

    if (!(isfile=getflag(argc, argv, "-f")))  isfile=0;
    if (!(verbose=getflag(argc, argv, "-v"))) verbose=0;
    if ((pstart=getflag(argc, argv, "-p")))   pstart++;
    else pstart = 1;

    if (!sgetpar(argc, argv, "-dir", &outdir)) outdir=NULL;

    /* Check parameters */
    if (fmin<=0.0 && fmax<=0.0) npoles = 0;
    if (wl<0.0) wl *= -1.0;
    if (!isfile && isatty(1)) \
        error("%s: can't write a SAC binary file to tty\n", argv[0]);


    /* Loop over SAC files */
    icomp = 0;
    nstat = 0;
    for (iarg=isfile+1; (iarg<argc && argv[iarg][0]!='-') || !isfile; iarg++) {

        /* Read SAC binary file */
        if (isfile) {
            if (!read_sacbin_file(argv[iarg], &hd, &data, &swap)) \
                error("%s: can't read %s\n", argv[0], argv[iarg]);
            if (verbose) \
                fprintf(stderr, "%s: station %d, file %d: %s\n", \
                    argv[0], nstat+1, icomp+1, argv[iarg]);
        }
        else {
            if (!read_sacbin(stdin, &hd, &data, &swap)) break;
        }

        /* Get required header values */
        dt = hd.delta;

        /* Check header values */
        if (dt<=0.0) {
            fprintf(stderr,"%s: invalid delta=%g, assuming delta=1.0\n", \
                argv[0], hd.delta);
            dt = 1.0;
        }


        /* Subtract mean value from trace data */
        if (hd.depmen!=SAC_HEADER_FLOAT_UNDEFINED) {
            for (it=0; it<hd.npts; it++) data[it] -= hd.depmen;
        }
        else {
            remove_mean(data, hd.npts);
        }        

        /* Butterworth low-cut and/or high-cut filter(s) */
        if (npoles>0) {
            butterworth_bandpass(npoles, fmin, fmax, dt, hd.npts, \
                zerophase, data);
        }


        /* Allocate space for three-component data */
        if (!icomp) {
            nt = hd.npts;
            if (!(data3 = fmalloc2(hd.npts, 3)) || \
                !(vdata = fmalloc2(9, hd.npts)) || \
                !(ddata = fmalloc2(3, hd.npts))) {
                error("%s: can't allocate space for 3-C data\n", argv[0]);
            }
        }

        /* Copy data */
        nt = (hd.npts>nt) ? nt : hd.npts;
        memcpy((void *)&hd3[icomp],  (const void *)&hd, sizeof(SACHEAD));
        memcpy((void *)data3[icomp], (const void *)data, nt*sizeof(float));

        /* Free space, and increment component counter */
        free(data);
        icomp++;

        /* Process three-component data */
        if (icomp==3) {
            icomp = 0;
            nstat++;

            /* Time window length in samples */
            iwl = NINT(wl/dt);
            if (iwl>nt-2) iwl = nt - 2;
            if (iwl<1)    iwl = 1;

            if (verbose) {
                fprintf(stderr,"%s: station %d: %d samples time window\n", \
                    argv[0], nstat, iwl);
            }

            /* P pick in samples (TO DO: implement field selection etc...) */
            if (hd3[0].a!=SAC_HEADER_FLOAT_UNDEFINED) {
                itp = NINT( (hd3[0].a - hd3[0].b) / dt);
            }
            else {
                itp = 0;
            }


            /* Covariance analysis */
            if (zeromean) do_eigen_zm(data3, nt, iwl, ddata, vdata);
            else          do_eigen(data3, nt, iwl, ddata, vdata);


            /* Allocate space for output trace */
            if (!(odata = fmalloc1(nt))) {
                error("%s: can't allocate space output trace\n", argv[0]);
            }


            /* Get adjusting parameter a and P polarization */
            a = ddata[iwl][0];
            if (itp>iwl && itp<nt-iwl) {
                itp = search_pvector1(vdata, ddata, nt, iwl, itp, a, &vp[0]);

                if (verbose) {
                    fprintf(stderr,"%s: station %d: ", argv[0], nstat);
                    fprintf(stderr, "a=%g theta=%g kp=%g\n", a, \
                        RAD2DEG*calc_theta(vp,0), calc_wuho_kp(ddata[itp], a));
                }
            }
            else {
                vp[0] = vp[1] = vp[2] = 0.0;
            }


            /* Loop over polarization attributes */
            for (j=pstart; (j<argc && argv[j][0]!='-'); j++) {
                isvalid = 1;

                if (STRCEQ(argv[j], "kp")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = calc_wuho_kp(ddata[it], a);
                    }
                }

                else if (STRCEQ(argv[j], "k2")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = calc_wuho_k2(ddata[it], a);
                    }
                }

                else if (STRCEQ(argv[j], "k1")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = calc_wuho_k1(vdata[it], vp);
                    }
                }

                else if (STRCEQ(argv[j], "k3")) {
                    odata[0] = 0.0;
                    for (it=1; it<iwl; it++) {
                        odata[it] = calc_wuho_k3(data3, vp, 0, it);
                    }
                    for (it=iwl; it<nt; it++) {
                        odata[it] = calc_wuho_k3(data3, vp, it-iwl, iwl);
                    }
                }

                else if (STRCEQ(argv[j], "ki")) {
                    for (it=0; it<nt; it++) {
                        odata[it]  = calc_wuho_k1(vdata[it], vp);
                        odata[it] *= calc_wuho_k2(ddata[it], a);
                        if (it<iwl)
                            odata[it] *= calc_wuho_k3(data3, vp, 0, it);
                        else
                            odata[it] *= calc_wuho_k3(data3, vp, it-iwl, iwl);

                        odata[it] *= odata[it];
                    }
                }

                /* Unknown polarization attribute name */
                else {
                    isvalid = 0;
                    if (nstat==1) {
                        fprintf(stderr, \
                            "\n%s: unknown attribute \"%s\", ignored\n",\
                            argv[0], argv[j]);
                    }
                }


                /* Set output SAC header */
                hd3[0].npts = nt;
                set_sac_depminmax(&hd3[0], odata);
                put_sac_hdr_string(&hd3[0], sac_hdr_index("kcmpnm"), argv[j]);

                /* Write SAC binary files... */
                if (isfile && isvalid) {
                    if (!out_file_name(argv[iarg-2], argv[j], outdir, &outfname)) \
                            error("%s: can't set output file name\n", argv[0]);

                    if (!write_sacbin_file(outfname, hd3[0], odata, swap))\
                        error("%s: can't write %s\n", argv[0], outfname);
                    free(outfname);
                }
                /* ... or write SAC binary to stdout */
                else if (isvalid) {
                    if (!write_sacbin(stdout, hd3[0], odata, swap)) \
                        error("%s: can't write SAC trace to stdout (%s)\n", \
                             argv[0], argv[j]);
                }
            }

            /* Free space of three-component data and output trace */
            free2((void *) data3);
            free2((void *) vdata);
            free2((void *) ddata);
            free((void *) odata);
        }
    }

    /* Diagnostic print */
    if (icomp) {
        fprintf(stderr, "%s: last %d trace%s skipped\n", \
            argv[0], icomp, (icomp==1)?"":"s");
    }

    return EXIT_SUCCESS;
}


/************************************************************************/
/* Functions used internally                                            */
/************************************************************************/

#define VECTORLENGTH(x,y,z)  ( sqrt((x)*(x) + (y)*(y) + (z)*(z)) )

float dot_product(float *v1, float *v2)
/************************************************************************
dot_product - compute dot product of two three-element vectors

Input:
v1      three-element vector (v1[3])
v2      three-element vector (v2[3])

Author: Nils Maercklin, 2009
*************************************************************************/
{
    register int i;
    float v=0.0;

    for (i=0; i<3; i++) v += v1[i]*v2[i];

    return v;
}



int search_pvector1(float **vdata, float **ddata, int nt, int iwl, \
    int itp, float a, float *vp)
/************************************************************************
search_pvector1 - search for P-wave vector around P-pick (EXPERIMENTAL)

Input:
vdata   two-dimensional array of principal eigenvectors (vdata[nt][3])
ddata   two-dimensional array of eigenvalues (ddata[nt][3])
nt      number of time samples in arrays vdata and ddata
iwl     analysis window length
itp     sample index of P-wave pick
a       adjusting parameter for attribute kp

Output:
vp      P-wave vector (vp[3]); function returns index of P-wave vector

Author: Nils Maercklin, 2009
*************************************************************************/
{
    register int it, j;
    int itpmax=itp;
    float valmax=0.0;
    float val=0.0;

    for (it=itp-iwl; it<itp+iwl && it<nt; it++) {
        if (it>=0) {
            val = fabs(vdata[it][0]) * calc_wuho_kp(ddata[it], a);
        }
        else {
            val = 0.0;
        }

        if (val>valmax) {
            valmax = val;
            itpmax = it;
            for (j=0; j<3; j++) vp[j] = vdata[it][j];
        }
    }

    return itpmax;
}



float calc_wuho_kp(float *d, float a)
/************************************************************************
calc_wuho_kp - compute Wu & Horiuchi polarization attribute kp

Input:
d       three-element array of eigenvalues (d[3])
a       adjusting parameter, largest eval of background cov. matrix

Reference:
Wu, C. and Horiuchi, S. (2008). Automatic determination of source 
    parameters of the 2007 Noto Hanto earthquake. Earth, Planets, and 
    Space, 60, 1053-1057. Equation 2.

Author: Nils Maercklin, 2009
*************************************************************************/
{
    float val;

    if (d[0]+a > 0.0) 
        val = (d[0] - d[1]) / (d[0] + a);
    else 
        val = 0.0;

    return val;
}



float calc_wuho_k1(float *v, float *vp)
/************************************************************************
calc_wuho_k1 - compute Wu & Horiuchi polarization attribute k1

Input:
v       three-element eigenvector V1 (v[3])
vp      normalized three-element P-wave vector (vp[3])

Reference:
Wu, C. and Horiuchi, S. (2008). Automatic determination of source 
    parameters of the 2007 Noto Hanto earthquake. Earth, Planets, and 
    Space, 60, 1053-1057. Equation 3.

Author: Nils Maercklin, 2009
*************************************************************************/
{
    float dp, val;

    dp = fabs(dot_product(v, vp));

    if      (dp<0.00001) val = 1.0;
    else if (dp>0.99999) val = 0.0;
    else val = (2.0/PI) * acos( dp );

    return val;
}



float calc_wuho_k2(float *d, float a)
/************************************************************************
calc_wuho_k2 - compute Wu & Horiuchi polarization attribute k2

Input:
d       three-element array of eigenvalues (d[3])
a       adjusting parameter, largest eval of background cov. matrix

Reference:
Wu, C. and Horiuchi, S. (2008). Automatic determination of source 
    parameters of the 2007 Noto Hanto earthquake. Earth, Planets, and 
    Space, 60, 1053-1057. Equation 4.

Author: Nils Maercklin, 2009
*************************************************************************/
{
    float nom, den, val;

    den = d[0]*d[0] + d[1]*d[1] + d[2]*d[2] + a*a;

    if (den) {
        nom = (d[0] - d[2])*(d[0] - d[2]) + (d[1] - d[2])*(d[1] - d[2]);
        val = nom / den;
    }
    else {
        val = 0.0;
    }

    return val;
}



float calc_wuho_k3(float **data3, float *vp, int it0, int iwl)
/************************************************************************
calc_wuho_k3 - compute Wu & Horiuchi polarization attribute k3

Input:
data3   three-component seismogram (data3[nt][3])
vp      normalized three-element P-wave vector (vp[3])
it0     first sample of analysis time window
iwl     length of analysis time window in seconds

Note:
This implementation is rather inefficient, because sums are built 
from scratch for each time sample.

Reference:
Wu, C. and Horiuchi, S. (2008). Automatic determination of source 
    parameters of the 2007 Noto Hanto earthquake. Earth, Planets, and 
    Space, 60, 1053-1057. Equation 5.

Author: Nils Maercklin, 2010
*************************************************************************/
{
    register int it;
    float u[3];
    float dp, uval, val;
    float num=0.0;
    float den=0.0;

    for (it=it0; it<it0+iwl; it++) {
        u[0] = data3[0][it];
        u[1] = data3[1][it];
        u[2] = data3[2][it];

        dp   = fabs(dot_product(u, vp));
        uval = VECTORLENGTH(u[0], u[1], u[2]);

        num += dp * dp;
        den += uval * uval;
    }

    if (iwl) {
        num /= (float)iwl;
        den /= (float)iwl;
    }


    if (iwl && den) {
        val = 1.0 - num/den;
    }
    else {
        val = 0.0;
    }

    return val;
}

/* END OF FILE */
