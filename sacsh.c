/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
SACSH - Set SAC header values, including different time format options
*************************************************************************

Author: Nils Maercklin,
    RISSC, University of Naples, Italy, October 2008

Version: 2008-10-31

Notes: Input files must be in SAC binary format (NVHDR=6). If multiple 
    traces are read from stdin, correct NPTS, LEVEN, and IFTYPE values
    are required. Only header values may be modified, not data sections.

Reference (SAC file format):
    SAC - Seismic Analysis Code, User's manual, 
    http://www.iris.edu/manuals/sac/manual.html (last checked: 2008-10-30)
*************************************************************************/

#include "nmsaclib.h"     /* SAC-related functions, and other stuff  */


/* Self-documentation */
const char *sdoc[] = {
"                                                                       ",
" SACSH - Set SAC Header field values                                   ",
"                                                                       ",
" Usage: sacsh [-k] field_name value [...] [options] -f sac_files       ",
"        sacsh <stdin [-k] field_name value [...] [options] > stdout    ",
"                                                                       ",
" Optional parameters:                                                  ",
"   -f        list of SAC binary files (alternatively read from stdin)  ",
"   -k        list of field names and values                            ",
"   -u        arbitrary string for header values to be undefined        ",
"   -d        flag: compute depmin, depmax, and depmen (if -f is set)   ",
"   -g        flag: recompute dist, gcarc, az, and baz                  ",
"   -v        flag: print warning messages                              ",
"   -y 1970   epoch year                                                ",
"                                                                       ",
" All SAC header fields and Z (zero time of file) can be set.           ",
" Times may be specified in seconds relative to Z (default), as epoch   ",
" time in seconds (prepend modifier +E), as the time of the day (+T),   ",
" or as a string of date and time (+H).                                 ",
" Symbolic values may be used for enum fields (modifier +N).            ",
"                                                                       ",
" Example:                                                              ",
" Set epicenter and origin time (as epoch time or as time string):      ",
"   sacsh evla 40.86 evlo 15.14 +eo 1223260046.2 -f *.sac               ",
"   sacsh evla 40.86 evlo 15.14 +ho 2008-10-06-02:27:26.2 -f *.sac      ",
"                                                                       ",
" NM, 2008-10-31                                                        ",
"                                                                       ",
NULL};


/* Definitions */
#define INDEX_ZTIME  1000

#ifndef WGS1984_A
#define WGS1984_A 6378137.0000
#endif
#ifndef WGS1984_F
#define WGS1984_F 0.003352810665
#endif

#ifdef IS_SACHDR_TIME
#undef IS_SACHDR_TIME
#endif
#define IS_SACHDR_TIME(id) (( ( (id)>=5 && (id)<=8 ) || \
    ( (id)>=SAC_HEADER_TMARK_POSITION && (id)<=SAC_HEADER_TMARK_POSITION+10 ))\
    || ( (id)==INDEX_ZTIME ) ? 1 : 0 )

#define IS_SACHDR_ENUM(id) \
    (( (id)>=SAC_HEADER_ENUM_MIN && (id)<=SAC_HEADER_ENUM_MAX ) ? 1 : 0 )

#define IS_OPTION(s)  (( STRCEQ(s,"-d") || STRCEQ(s,"-f") || \
    STRCEQ(s,"-g") || STRCEQ(s,"-k") || STRCEQ(s,"-u") || STRCEQ(s,"-v") || \
    STRCEQ(s,"-y")) ? 1 : 0)


/* Prototypes of functions used internally */
int set_time(SACHEAD *hd, int index, char modif, int eyear, char *undef, \
    char *val);
double decode_time(char *val, char modif, int eyear);
int check_zdate(SACHEAD hd);
int check_ztime(SACHEAD hd);
int comp_azim_dist(SACHEAD *hd);
int comp_minmax(FILE *fp, SACHEAD *hd, int swap);


int 
main(int argc, char **argv)
{
    /* Variables */
    FILE *fp=NULL;      /* input file pointer */
    register int i,j,iarg;  /* loop indices */
    SACHEAD hd;         /* SAC header */
    int isfile;         /* flag: file names (position of "-f" option) */
    int kstart=1;       /* flag: key (position of "-k" option) */
    int verbose=0;      /* flag: print warnings */
    int issac=0;        /* flag: SAC header appears to be valid */
    int swap=0;         /* flag: header byte-swapping */
    int lcalda=0;       /* flag: recompute gcarc, az, baz, dist */
    int minmax=0;       /* flag: compute depmin, depmax, depmen (if "-f") */
    char *undef=NULL;   /* output string for undefined values */
    char *key=NULL;     /* SAC header field name */
    int index;          /* SAC header field index */
    int npts;           /* number of samples in trace */
    int ival;           /* value of integer field */
    float fval;         /* value of float field */
    int ntr=0;          /* input trace counter */
    char modif=' ';     /* format modifier */
    int eyear=1970;     /* epoch year */

    /* Print documentation */
    if (argc==1) {
        i=0;
        while (sdoc[i]) fprintf(stderr, "%s\n", sdoc[i++]);
        return (EXIT_FAILURE);
    }

    /* Command line options */
    if (!(igetpar(argc, argv, "-y", &eyear))) eyear=1970; 
    if (!(sgetpar(argc, argv, "-u", &undef))) undef=SAC_HEADER_UNDEFINED;
    if (!(isfile=getflag(argc, argv, "-f")))  isfile=0;
    if (!(verbose=getflag(argc, argv, "-v"))) verbose=0;
    if (!(lcalda=getflag(argc, argv, "-g")))  lcalda=0;
    if (!(minmax=getflag(argc, argv, "-d")))  minmax=0;
    if ((kstart=getflag(argc, argv, "-k")))   kstart++;
    else kstart = 1;

    /* Validation */
    if (minmax && !isfile) {
        minmax = 0;
        fprintf(stderr, "%s: option \"-d\" ignored for stdin\n", argv[0]);
    }
    if (strlen(undef)==0) undef=NULL;


    /* Loop over SAC files (or read from stdin) */
    for (iarg=isfile+1; (iarg<argc && !IS_OPTION(argv[iarg])) || !isfile; \
         iarg++) {

        /* Read SAC binary file header */
        if (isfile) {
            if ((fp=fopen(argv[iarg], "r+")) && \
                (fread(&hd, sizeof(SACHEAD), 1, fp) == 1)) {
                issac = 1;
            }
            else {
                issac = 0;
                if (verbose) {
                    fprintf(stderr, "%s: can't read %s\n", argv[0], argv[iarg]);
                }
            }
        }
        else {
            if (fread(&hd, sizeof(SACHEAD), 1, stdin) != 1) break;
            issac = 1;
        }

        /* Byte-swapping, and check, if hd appears to be SAC header */
        swap = 0;
        if (issac && NEED_BYTESWAP(hd.nvhdr)) {
            swab4((char*) &hd, SAC_HEADER_SIZE_NUMBERS);
            swap = 1;

            if (NEED_BYTESWAP(hd.nvhdr)) {
                issac = 0;
                if (verbose) {
                    fprintf(stderr, "%s: invalid header in %s\n", \
                        argv[0], (isfile) ? argv[iarg] : "stdin" );
                }
            }
        }

        /* Increment trace counter */
        if (issac) ntr++;


        /* Loop over specified SAC header fields */
        for (j=kstart; (j<argc-1 && !IS_OPTION(argv[j]) && issac); j+=2) {

            /* Get modifier and SAC header field name */
            if (argv[j][0]=='+' && strlen(argv[j])>2) {
                modif = argv[j][1];
                key   = argv[j]+2;
            }
            else {
                modif = ' ';
                key   = argv[j];
            }

            /* Get index of header field (negative, if invalid) */
            if (STRCEQ(key,"z")) {
                index = INDEX_ZTIME;
            }
            else {
                index = sac_hdr_index(key);
            }

            /* Ignore invalid field names */
            if (index < 0) {
                if (verbose && ntr==1) {
                    fprintf(stderr, "%s: invalid field name %s\n", \
                        argv[0], key);
                }
            }

            /* Process time fields (special input format otions) */
            else if (IS_SACHDR_TIME(index)) {
                if (!set_time(&hd, index, modif, eyear, undef, argv[j+1]) \
                    && verbose) {
                    fprintf(stderr, "%s: can't set %s in %s\n", \
                        argv[0], key, (isfile) ? argv[iarg] : "stdin");
                }
            }

            /* Process other float fields */
            else if (IS_SACHDR_FLOAT(index)) {
                if (STRCEQ(argv[j+1], undef)) {
                    fval = SAC_HEADER_FLOAT_UNDEFINED;
                }
                else {
                    fval = atof(argv[j+1]);
                }
                put_sac_hdr_float(&hd, index, fval);
            }

            /* Process integer fields, including enums */
            else if (IS_SACHDR_INT(index)) {
                if (STRCEQ(argv[j+1], undef)) {
                    ival = SAC_HEADER_INT_UNDEFINED;
                }
                else if (IS_SACHDR_ENUM(index) && (modif=='n' || modif=='N')) {
                    ival = sac_hdr_enum(argv[j+1]);
                    if (ival<0) {
                        if (verbose && ntr==1) {
                            fprintf(stderr, "%s: invalid enum name %s\n", \
                                argv[0], argv[j+1]);
                        }
                        ival = get_sac_hdr_int(&hd, index);
                    }
                }
                else {
                    ival = atoi(argv[j+1]);
                }
                put_sac_hdr_int(&hd, index, ival);
            }

            /* Process character fields */
            else {
                if (STRCEQ(argv[j+1], undef)) {
                    put_sac_hdr_string(&hd, index, SAC_HEADER_CHAR_UNDEFINED);
                }
                else {
                    put_sac_hdr_string(&hd, index, argv[j+1]);
                }
            }
        } /* end of key-loop */


        /* Recompute azimuth and distance header values */
        if (lcalda) {
            if (!comp_azim_dist(&hd) && verbose) {
                fprintf(stderr, "%s: undefined coordinates in %s\n", \
                    argv[0], (isfile) ? argv[iarg] : "stdin");
            }
        }

        /* Compute depmin, depmax, and depmen */
        if (minmax) {
            if (!comp_minmax(fp, &hd, swap) && verbose) {
                fprintf(stderr, "%s: can't compute min/max/mean for %s\n", \
                    argv[0], (isfile) ? argv[iarg] : "stdin");
            }
        }


        /* Remember npts and swap header bytes, if necessary */
        npts = hd.npts;
        if (swap) swab4((char*) &hd, SAC_HEADER_SIZE_NUMBERS);


        /* Write modified header to specified file, and close file */
        if (isfile && issac) {
            rewind(fp);
            if (fwrite(&hd, sizeof(SACHEAD), 1, fp) != 1) {
                fprintf(stderr, "%s: write failed for %s\n", \
                    argv[0], argv[iarg]);
            }
            fclose(fp);
        }

        /* Write SAC trace to stdout (if input from stdin) */
        else {
            /* SAC header */
            if (fwrite(&hd, sizeof(SACHEAD), 1, stdout) != 1) {
                error("%s: write error on stdout\n", argv[0]);
            }

            /* Trace data (1 or 2 sections) */
            ival = (!hd.leven || hd.iftype==IRLIM || hd.iftype==IAMPH) ? 2 : 1;
            for (j=0; j<ival*npts; j++) {
                if (fread(&fval, sizeof(float), 1, stdin)  != 1) \
                    error("%s: read error from stdin\n", argv[0]);
                if (fwrite(&fval, sizeof(float), 1, stdout) != 1) \
                    error("%s: write error on stdout\n", argv[0]);
            }
        }

    } /* end of file-loop */


    return EXIT_SUCCESS;
}


/************************************************************************/
/* Functions used internally                                            */
/************************************************************************/

int check_zdate(SACHEAD hd)
/************************************************************************
check_zdate - check definition of file reference date (zero time)

Input:
hd      SAC header structure   (see SACHEAD typedef)

Output: returns 1 (true) or 0 (false)

Author: Nils Maercklin, 2008-10-05
*************************************************************************/
{
    if (hd.nzyear != SAC_HEADER_INT_UNDEFINED && \
        hd.nzjday != SAC_HEADER_INT_UNDEFINED) {
        return 1;
    }
    else {
        return 0;
    }
}



int check_ztime(SACHEAD hd)
/************************************************************************
check_ztime - check definition of file reference time (zero time)

Input:
hd      SAC header structure   (see SACHEAD typedef)

Output: returns 1 (true) or 0 (false)

Author: Nils Maercklin, 2008-10-05
*************************************************************************/
{
    if (hd.nzhour != SAC_HEADER_INT_UNDEFINED && \
        hd.nzmin  != SAC_HEADER_INT_UNDEFINED && \
        hd.nzsec  != SAC_HEADER_INT_UNDEFINED && \
        hd.nzmsec != SAC_HEADER_INT_UNDEFINED) {
        return 1;
    }
    else {
        return 0;
    }
}



int set_time(SACHEAD *hd, int index, char modif, int eyear, char *undef, \
    char *val)
/************************************************************************
set_time - set SAC time header field value

Input:
hd      pointer to SAC header structure (see SACHEAD typedef)
index   index of header field  (time), or INDEX_ZTIME
modif   output format modifier (see switch statement)
eyear   epoch year
undef   string used for undefined header fields
val     new value of time field

Output: modified SAC header structure (function returns 0 on failure)

Author: Nils Maercklin, 2008-10-31
*************************************************************************/
{
    double etime=0.0;
    double eztime=0.0;
    float sec=0.0;    
    int day=0;
    int year=0;
    int hour=0;
    int minute=0;


    /* Set file reference time (zero time) */
    if (index == INDEX_ZTIME && STRCEQ(val,undef)) {
        hd->nzyear = SAC_HEADER_INT_UNDEFINED;
        hd->nzjday = SAC_HEADER_INT_UNDEFINED;
        hd->nzhour = SAC_HEADER_INT_UNDEFINED;
        hd->nzmin  = SAC_HEADER_INT_UNDEFINED;
        hd->nzsec  = SAC_HEADER_INT_UNDEFINED;
        hd->nzmsec = SAC_HEADER_INT_UNDEFINED;    
    }
    else if (index == INDEX_ZTIME) {
        switch (modif) {
            case 'H':    /* Human date and time */
            case 'h':
                etime = decode_time(val, modif, eyear);
                etime_to_date(etime, eyear, &year, &day, &hour, &minute, &sec);
                break;
            case 'T':    /* Time of day (hour, minute, second) */
            case 't':
                etime = decode_time(val, modif, eyear);
                etime_to_date(etime, eyear, &year, &day, &hour, &minute, &sec);
                year  = hd->nzyear;
                day   = hd->nzjday;
                break;
            case 'E':    /* Epoch time in seconds */
            case 'e':
                etime = atof(val);
                break;
            default:     /* Time shift by val seconds */
                if (!check_zdate(*hd) || !check_ztime(*hd)) {
                    return 0;
                }
                else {
                    sec = ((float) hd->nzsec) + 0.001*((float) hd->nzmsec) + \
                        atof(val);
                    etime = date_to_etime(eyear, hd->nzyear, hd->nzjday, \
                        hd->nzhour, hd->nzmin, sec);
                }
                etime_to_date(etime, eyear, &year, &day, &hour, &minute, &sec);
                break;
        }

        hd->nzyear = year;
        hd->nzjday = day;
        hd->nzhour = hour;
        hd->nzmin  = minute;
        hd->nzsec  = (int) floor(sec);
        hd->nzmsec = NINT(1000.0*(sec-floor(sec)));
    }

    /* Set other time fields */
    else if (STRCEQ(val,undef)) {
        put_sac_hdr_float(hd, index, SAC_HEADER_FLOAT_UNDEFINED);
    }
    else {
        switch (modif) {
            case 'H':    /* Human date and time */
            case 'h':
                etime = decode_time(val, modif, eyear);
                if (!check_zdate(*hd) || !check_ztime(*hd)) {
                    put_sac_hdr_float(hd, index, SAC_HEADER_FLOAT_UNDEFINED);
                    return 0;
                }
                else {
                    sec    = ((float)hd->nzsec) + 0.001*((float)hd->nzmsec);
                    eztime = date_to_etime(eyear, hd->nzyear, hd->nzjday, \
                        hd->nzhour, hd->nzmin, sec);

                    put_sac_hdr_float(hd, index, (float)(etime-eztime));
                }
                break;
            case 'T':    /* Time of day (hour, minute, second) */
            case 't':
                etime = decode_time(val, modif, eyear);
                if (!check_ztime(*hd)) {
                    put_sac_hdr_float(hd, index, SAC_HEADER_FLOAT_UNDEFINED);
                    return 0;
                }
                else {
                    sec    = ((float)hd->nzsec) + 0.001*((float)hd->nzmsec);
                    eztime = date_to_etime(eyear, eyear, 1, hd->nzhour, \
                        hd->nzmin, sec);

                    put_sac_hdr_float(hd, index, (float)(etime-eztime));
                }
                break;
            case 'E':    /* Epoch time in seconds */
            case 'e':
                if (!check_zdate(*hd) || !check_ztime(*hd)) {
                    put_sac_hdr_float(hd, index, SAC_HEADER_FLOAT_UNDEFINED);
                    return 0;
                }
                else {
                    sec    = ((float)hd->nzsec) + 0.001*((float)hd->nzmsec);
                    eztime = date_to_etime(eyear, hd->nzyear, hd->nzjday, \
                        hd->nzhour, hd->nzmin, sec);
                    etime  = atof(val);

                    put_sac_hdr_float(hd, index, (float)(etime-eztime));
                }
                break;
            default:     /* Seconds relative to reference time */
                put_sac_hdr_float(hd, index, atof(val));
                break;
        }
    }

    return 1;
}



double decode_time(char *val, char modif, int eyear)
/************************************************************************
decode_time - decode a time string and return time in seconds

Input:
val     time string
modif   input format modifier (see switch statement)
eyear   epoch year

Output:
etime   epoch time in seconds, or seconds of the day
        (function returns 0 on failure, i.e. if month is invalid)

Author: Nils Maercklin, 2008-10-31
*************************************************************************/
{
    double etime=0.0;
    char *ptr=NULL;
    int year=eyear;
    int month=1;
    int day=1;
    int hour=0;
    int minute=0;
    float sec=0.0;

    ptr = val;
    switch (modif) {
        case 'H':    /* Human date and time */
        case 'h':
            if (strlen(ptr)) year   = strtol(ptr,   &ptr, 10);
            if (strlen(ptr)) month  = strtol(ptr+1, &ptr, 10);
            while (month<=0) {
                year--;
                month+=12;
            }
            if (month>12) {
                year += month/12;
                month = month%12;
            }
            if (strlen(ptr)) day    = strtol(ptr+1, &ptr, 10);
            if (strlen(ptr)) hour   = strtol(ptr+1, &ptr, 10);
            if (strlen(ptr)) minute = strtol(ptr+1, &ptr, 10);
            if (strlen(ptr)) sec    = strtod(ptr+1, &ptr);

            day   = hdate_to_jday(year, month, day);
            etime = date_to_etime(eyear, year, day, hour, minute, sec);
            break;
        case 'T':    /* Time of day (hour, minute, second) */
        case 't':
            if (strlen(ptr)) hour   = strtol(ptr,   &ptr, 10);
            if (strlen(ptr)) minute = strtol(ptr+1, &ptr, 10);
            if (strlen(ptr)) sec    = strtod(ptr+1, &ptr);

            etime = date_to_etime(eyear, eyear, 1, hour, minute, sec);
            break;
        case 'E':    /* Epoch time in seconds */
        case 'e':
        default:
            etime = atof(val);
            break;
    }

    return etime;
}



int comp_azim_dist(SACHEAD *hd)
/************************************************************************
comp_azim_dist - recompute azimuth and distance header values

Input:
hd      pointer to SAC header structure (see SACHEAD typedef)

Output: modified gcarc, dist, az, baz in SAC header structure
        (function returns 0 on failure)

Author: Nils Maercklin, 2008-10-31
*************************************************************************/
{
    double evla, evlo, stla, stlo;
    double gcarc, dist, az, baz;

    if (hd->evla != SAC_HEADER_INT_UNDEFINED && \
        hd->evlo != SAC_HEADER_INT_UNDEFINED && \
        hd->stla != SAC_HEADER_INT_UNDEFINED && \
        hd->stlo != SAC_HEADER_INT_UNDEFINED) {

        evla = (double) hd->evla;
        evlo = (double) hd->evlo;
        stla = (double) hd->stla;
        stlo = (double) hd->stlo;

        gcarc = gc_delta(evla, evlo, stla, stlo);
        dist  = 0.001 * gc_dist(evla, evlo, stla, stlo, WGS1984_A, WGS1984_F);
        az    = gc_azimuth(evla, evlo, stla, stlo);
        baz   = gc_azimuth(stla, stlo, evla, evlo);

        hd->gcarc = (float) gcarc;
        hd->dist  = (float) dist;
        hd->az    = (float) az;
        hd->baz   = (float) baz;

        return 1;
    }
    else {
        return 0;
    }

    return 0;
}



int comp_minmax(FILE *fp, SACHEAD *hd, int swap)
/************************************************************************
comp_minmax - compute minimum, maximum, and mean of dependent variable

Input:
fp      pointer to SAC file (should be a disk file)
hd      pointer to SAC header structure (see SACHEAD typedef)
swap    byte-swapping flag (1 = swap data bytes, 0 = do not)

Output: modified depmin, depmax, depmen in SAC header structure
        (function returns 0 on failure, i.e. a read error)

Author: Nils Maercklin, 2008-10-31
*************************************************************************/
{
    FILE *tmpfp=NULL;
    double dmin=0.0;
    double dmax=0.0;
    double dmean=0.0;
    float  fval=0.0;
    int    npts=0;
    register int i;

    tmpfp = fp;
    npts  = hd->npts;

    if (npts==SAC_HEADER_INT_UNDEFINED) {
        return 0;
    }
    else {
        for (i=0; i<npts; i++) {
            if (fread(&fval, sizeof(float), 1, tmpfp) != 1) return 0;
            if (swap) swab4((char*) &fval, sizeof(float));

            dmean += (double) fval;

            if (!i || (double)fval > dmax) dmax = (double) fval;
            if (!i || (double)fval < dmin) dmin = (double) fval;
        }
        dmean /= (double) npts;

        hd->depmin = (float) dmin;
        hd->depmax = (float) dmax;
        hd->depmen = (float) dmean;
    }

    return 1;
}

/* END OF FILE */
