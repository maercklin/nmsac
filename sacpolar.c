/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
SACPOLAR - SAC POLARization analysis of three-component (3-C) data
*************************************************************************

Author: Nils Maercklin,
    RISSC, University of Naples, Italy, September 2009

Version: 2011-06-12

Modifications:
  2009-09-11 (NM): First version of this code
  2010-01-27 (NM): Added output of attribute name, header KCMPNM
  2010-02-04 (NM): Added normalized inclination angles inc1, inc3
  2010-03-17 (NM): Added option -z for fast eigenanalysis of zero-mean 
                   data, and enabled optional zerophase bandpass
  2011-06-12 (NM): Added user-specified output directory

Notes: Input files must be evenly-sampled time series files in SAC 
    binary format (assuming NVHDR=6, IFTYPE=ITIME (1), LEVEN=1).
    Three subsequent files/traces are considered as one 3-C dataset with 
    common trace length, start time, and sampling rate.

    No station/component consistency checks are made (yet), and the 
    code may contain some (undocumented) experimental features.

    Polarization attribute definitions are given in Maercklin (1999) and 
    in my SUPOLAR manual available at http://purl.org/net/nils/man/supolar.


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
Samson, J. C. (1973). Descriptions of the polarization states of vector 
    processes: applications to ULF magnetic fields. 
    Geophys. J. R. Astr. Soc., 34(4), 403-419.
*************************************************************************/

#include "nmsaclib.h"     /* SAC-related functions, and other stuff */


/* Self-documentation */
const char *sdoc[] = {
"                                                                       ",
" SACPOLAR - SAC Polarization analysis of three-component data          ",
"                                                                       ",
" Usage: sacpolar [-p] attribute [...] [parameters] -f sac_files        ",
"        sacpolar <stdin [-p] attribute [...] [parameters] > stdout     ",
"                                                                       ",
" Parameters and defaults:                                              ",
"   -f        list of SAC binary waveform files (or read from stdin)    ",
"   -p        list of one or more of the pol. attributes listed below   ",
"   -w   0.5  correlation time window length in seconds                 ",
"   -q   1.0  contrast parameter of rectilinearity RL                   ",
"   -z        flag: assume zero mean in correlation windows  (faster)   ",
"   -v        flag: verbose operation                                   ",
"   -dir      output directory (default is same as input directory)     ",
"                                                                       ",
"   -b1    0  Butterworth low-cut frequency in Hz  (0 = no filter)      ",
"   -b2    0  Butterworth high-cut frequency in Hz (0 = no filter)      ",
"   -bp    3  number of poles of Butterworth filters                    ",
"   -bz       flag: apply zerophase (two-pass) Butterworth filter       ",
"                                                                       ",
" Polarization attributes:                                              ",
"   rl, rl2   rectilinearity RL                             [0, 1]      ",
"   tau       global polarization parameter tau             [0, 1]      ",
"   l1        linearity coefficient l1                      [0, 1]      ",
"   f1, pln   flatness coefficient f1 and planarity         [0, 1]      ",
"   inc1,inc3 normalized long- and short-axis inclination   [0, 1]      ",
"   theta     vertical polarization angle (incidence)  [   0,  90] deg  ",
"   phi, phi1 horizontal polarization angle            [ -90,  90] deg  ",
"   phi2      horizontal polarization angle            [-180, 180] deg  ",
"   phi3      horizontal polarization angle            [   0, 360] deg  ",
"   er        eigenresultant (polarization amplitude)  [   0, inf]      ",
"                                                                       ",
" The polarization attributes are computed in a moving time window from ",
" the eigenvalues and the principal eigenvector of the three-component  ",
" covariance matrix.                                                    ",
" Three subsequent files/traces are considered as one three-component   ",
" dataset with common trace length, start time, and sampling rate.      ",
" For correct angles, the vertical component has to be the first trace. ",
" Each output file has a suffix corresponding to the attribute name.    ",
"                                                                       ",
" NM, 2011-06-12                                                        ",
"                                                                       ",
NULL};


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

    float wl=0.0;        /* polarization time window length in seconds */
    int iwl=0;           /* polarization time window length in samples */

    int pstart=1;        /* start of polarization attribute list */
    float rlq=0.0;       /* contrast factor of rectilinearity RL */


    /* Print documentation */
    if (argc==1) {
        i=0;
        while (sdoc[i]) fprintf(stderr, "%s\n", sdoc[i++]);
        return (EXIT_FAILURE);
    }

    /* Command line parameters */
    if (!fgetpar(argc, argv, "-w",  &wl))      wl=0.5;
    if (!fgetpar(argc, argv, "-q",  &rlq))     rlq=1.0;

    if (!fgetpar(argc, argv, "-b1", &fmin))    fmin=0.0;
    if (!fgetpar(argc, argv, "-b2", &fmax))    fmax=0.0;
    if (!igetpar(argc, argv, "-bp", &npoles))  npoles=3;
    if (!(zerophase=getflag(argc, argv, "-bz"))) zerophase=0;

    if (!(zeromean=getflag(argc, argv, "-z")) && \
        !(zeromean=getflag(argc, argv, "-zm")))  zeromean=0;

    if (!(isfile=getflag(argc, argv, "-f")))   isfile=0;
    if (!(verbose=getflag(argc, argv, "-v")))  verbose=0;
    if ((pstart=getflag(argc, argv, "-p")))    pstart++;
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


            /* Covariance analysis */
            if (zeromean) do_eigen_zm(data3, nt, iwl, ddata, vdata);
            else          do_eigen(data3, nt, iwl, ddata, vdata);


            /* Allocate space for output trace */
            if (!(odata = fmalloc1(nt))) {
                error("%s: can't allocate space for output trace\n", argv[0]);
            }

            /* Loop over polarization attributes */
            for (j=pstart; (j<argc && argv[j][0]!='-'); j++) {
                isvalid = 1;

                /* Rectilinearity RL */
                if (STRCEQ(argv[j], "rl") || STRCEQ(argv[j], "rl1")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = calc_rl(ddata[it], rlq, 0);
                    }
                }
                else if (STRCEQ(argv[j], "rl2")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = calc_rl(ddata[it], rlq, 1);
                    }
                }

                /* Planarity (Jurkevics, 1988) */
                else if (STRCEQ(argv[j], "pln")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = calc_plan(ddata[it]);
                    }
                }

                /* Flatness coefficient f1 (Benhama et al., 1988) */
                else if (STRCEQ(argv[j], "f1")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = calc_f1(ddata[it]);
                    }
                }

                /* Linearity coefficient l1  */
                else if (STRCEQ(argv[j], "l1")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = calc_l1(ddata[it]);
                    }
                }

                /* Global polarization parameter tau (Samson, 1973)  */
                else if (STRCEQ(argv[j], "tau")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = calc_tau(ddata[it]);
                    }
                }

                /* Ellipticities */
                else if (STRCEQ(argv[j], "e21") || \
                         STRCEQ(argv[j], "e31") || \
                         STRCEQ(argv[j], "e32")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = calc_ellip(ddata[it], \
                            argv[j][1]-'1', argv[j][2]-'1');
                    }
                }

                /* Eigenresultant (amplitude in polarization direction) */
                else if (STRCEQ(argv[j], "er")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = calc_er(ddata[it]);
                    }
                }

                /* Normalized long- and short-axis inclination angles */
                else if (STRCEQ(argv[j], "inc1") || \
                         STRCEQ(argv[j], "inc2") || \
                         STRCEQ(argv[j], "inc3") ) {
                    int ic = atoi(argv[j]+3) - 1;
                    for (it=0; it<nt; it++) {
                        odata[it] = calc_norminc(vdata[it]+(3*ic));
                    }
                }

                /* Incident angle thata */
                else if (STRCEQ(argv[j], "theta") || STRCEQ(argv[j], "inc")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = RAD2DEG * calc_theta(vdata[it], 0);
                    }
                }

                /* Azimuth phi */
                else if (STRCEQ(argv[j], "phi") || STRCEQ(argv[j], "phi1")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = RAD2DEG * calc_phi(vdata[it], 1);
                    }
                }
                else if (STRCEQ(argv[j], "phi2")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = RAD2DEG * calc_phi(vdata[it], 2);
                    }
                }
                else if (STRCEQ(argv[j], "phi3")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = RAD2DEG * calc_phi(vdata[it], 3);
                    }
                }

                /* Eigenvalues and eigenvectors (undocumented feature) */
                else if (STRCEQ(argv[j], "d1") || \
                         STRCEQ(argv[j], "d2") || \
                         STRCEQ(argv[j], "d3")) {
                    for (it=0; it<nt; it++) {
                        odata[it] = ddata[it][argv[j][1]-'1'];
                    }
                }
                else if (STRCEQ(argv[j], "v11") || STRCEQ(argv[j], "v21") || \
                         STRCEQ(argv[j], "v31") || STRCEQ(argv[j], "v12") || \
                         STRCEQ(argv[j], "v22") || STRCEQ(argv[j], "v32") || \
                         STRCEQ(argv[j], "v13") || STRCEQ(argv[j], "v23") || \
                         STRCEQ(argv[j], "v33")) {
                    int ic = argv[j][2] - '1' + 3*(argv[j][1] - '1');
                    for (it=0; it<nt; it++) {
                        odata[it] = vdata[it][ic];
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

/* END OF FILE */
