/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
SACLH - List SAC header values, including different time format options
*************************************************************************

Author: Nils Maercklin,
    RISSC, University of Naples, Italy, October 2008

Version: 2008-11-05

Notes: Input files must be in SAC binary format (NVHDR=6). If multiple 
    traces are read from stdin, correct NPTS, LEVEN, and IFTYPE values
    are required.

Reference (SAC file format):
    SAC - Seismic Analysis Code, User's manual, 
    http://www.iris.edu/manuals/sac/manual.html (last checked: 2008-10-30)
*************************************************************************/

#include "nmsaclib.h"     /* SAC-related functions, and other stuff */


/* Self-documentation */
const char *sdoc[] = {
"                                                                       ",
" SACLH - List SAC Header field values                                  ",
"                                                                       ",
" Usage: saclh [-k] field_names [options] -f sac_files                  ",
"        saclh <stdin [-k] field_names [options]                        ",
"                                                                       ",
" Optional parameters:                                                  ",
"   -f        list of SAC binary files (alternatively read from stdin)  ",
"   -k        list of SAC header field names to be printed to stdout    ",
"   -s \\t     output field separator (default: tab-separated list)     ",
"   -b 1      remove some blanks from character fields                  ",
"             (0: none, 1: at end, 2: at beginning and end)             ",
"   -u        arbitrary output string for undefined header values       ",
"             (default: UNDEF, \"+\" prints the true header value)      ",
"   -h        flag: include column description line in output list      ",
"   -v        flag: print warning messages                              ",
"   -y 1970   epoch year                                                ",
"                                                                       ",
" All SAC header fields, and NTR (trace number), FILE (SAC file name),  ",
" and Z (zero time of file) can be specified.                           ",
" Times are printed in seconds relative to Z (default), as epoch time   ",
" in seconds (prepend modifier +E), or as a human-readable string of    ",
" date and time (+H, +N) or of the time of the day only (+T).           ",
" For enum fields symbolic values (+N) or a description (+D) may be     ",
" printed instead of the corresponding integer value.                   ",
"                                                                       ",
" Example:                                                              ",
" List station names and first-arrivals (header value and time string): ",
"   saclh kstnm a +ha -f *.sac                                          ",
"                                                                       ",
" NM, 2008-11-05                                                        ",
"                                                                       ",
NULL};


/* Definitions */
#define INDEX_ZTIME  1000

#ifdef IS_SACHDR_TIME
#undef IS_SACHDR_TIME
#endif
#define IS_SACHDR_TIME(id) (( ( (id)>=5 && (id)<=8 ) || \
    ( (id)>=SAC_HEADER_TMARK_POSITION && (id)<=SAC_HEADER_TMARK_POSITION+10 ))\
    || ( (id)==INDEX_ZTIME ) ? 1 : 0 )

#define IS_SACHDR_ENUM(id) \
    (( (id)>=SAC_HEADER_ENUM_MIN && (id)<=SAC_HEADER_ENUM_MAX ) ? 1 : 0 )


/* Prototypes of functions used internally */
void print_enum(SACHEAD hd, int index, char modif, char *undef);
void print_time(SACHEAD hd, int index, char modif, int eyear, char *undef);
int check_zdate(SACHEAD hd);
int check_ztime(SACHEAD hd);

int 
main(int argc, char **argv)
{
    /* Variables */
    FILE *fp=NULL;      /* input file pointer */
    register int i,j,k,iarg; /* loop indices */
    SACHEAD hd;         /* SAC header */
    int isfile;         /* flag: file names (position of "-f" option) */
    int headerline=0;   /* flag: output header line */
    int kstart=1;       /* flag: key (position of "-k" option) */
    int verbose=0;      /* flag: print warnings */
    int issac=0;        /* flag: SAC header appears to be valid */
    int rmblank=2;      /* remove-whitespace option for character fields */
    char *fieldsep=NULL;/* output field separator */
    char *undef=NULL;   /* output string for undefined values */
    char *key=NULL;     /* SAC header field name */
    int index;          /* SAC header field index */
    int ival;           /* value of integer field */
    float fval;         /* value of float field */
    char string[17];    /* string buffer */
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
    if (!(igetpar(argc, argv, "-y", &eyear)))    eyear=1970; 
    if (!(igetpar(argc, argv, "-b", &rmblank)))  rmblank=1; 
    if (!(sgetpar(argc, argv, "-s", &fieldsep))) fieldsep="\t"; 
    if (!(sgetpar(argc, argv, "-u", &undef)))    undef=SAC_HEADER_UNDEFINED;
    if (strlen(undef)==0 || (strlen(undef)>0 && undef[0]=='+')) undef=NULL;

    if (!(isfile=getflag(argc, argv, "-f")))     isfile=0;
    if (!(headerline=getflag(argc, argv, "-h"))) headerline=0;
    if (!(verbose=getflag(argc, argv, "-v")))    verbose=0;
    if ((kstart=getflag(argc, argv, "-k")))      kstart++;
    else kstart = 1;


    /* Print field names as header line */
    if (headerline) {
        for (j=kstart,k=0; (j<argc && argv[j][0]!='-'); j++,k++) {
            if (argv[j][0]=='+' && strlen(argv[j])>2) key = argv[j]+2;
            else key = argv[j];
            printf("%s%d=%s", (k) ? fieldsep : "", k+1, key);
        }
        printf("\n");
     }


    /* Loop over SAC files (or read from stdin) */
    for (iarg=isfile+1; (iarg<argc && argv[iarg][0]!='-') || !isfile; iarg++) {

        /* Read SAC binary file header */
        if (isfile) {
            if ((fp=fopen(argv[iarg], "r")) && \
                (fread(&hd, sizeof(SACHEAD), 1, fp) == 1)) {
                issac = 1;
                fclose(fp);
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
        if (issac && NEED_BYTESWAP(hd.nvhdr)) {
            swab4((char*) &hd, SAC_HEADER_SIZE_NUMBERS);

            if (NEED_BYTESWAP(hd.nvhdr)) {
                issac = 0;
                if (verbose) {
                    fprintf(stderr, "%s: invalid header in %s\n", \
                        argv[0], (isfile) ? argv[iarg] : "stdin" );
                }
            }
        }


        /* Increment trace counter and print newline */
        if (issac && ++ntr>1) printf("\n");


        /* Loop over specified SAC header fields */
        for (j=kstart,k=0; (j<argc && argv[j][0]!='-' && issac); j++,k++) {

            /* Get modifier and SAC header field name */
            if (argv[j][0]=='+' && strlen(argv[j])>2) {
                modif = argv[j][1];
                key   = argv[j]+2;
            }
            else {
                modif = ' ';
                key   = argv[j];
            }

            /* Get index of header field */
            index = sac_hdr_index(key);

            /* Print field separator */
            if (k) printf("%s", fieldsep);


            /* Process "pseudo fields" and invalid field names */
            if (index < 0) {
                if (STRCEQ(key,"ntr")) {
                    printf("%d", ntr);
                }
                else if (STRCEQ(key,"file")) {
                    printf("%s", (isfile) ? argv[iarg] : "stdin");
                }
                else if (STRCEQ(key,"z")) {
                    print_time(hd, INDEX_ZTIME, modif, eyear, undef);
                }
                else if (STRCEQ(key,"kztime")) {
                    print_time(hd, INDEX_ZTIME, 't', eyear, undef);
                }
                else {
                    printf("%s", "(NULL)");
                    if (ntr==1 && verbose) {
                        fprintf(stderr, "%s: invalid field name %s\n", \
                            argv[0], key);
                    }
                }
            }

            /* Process time and enum fields (special formatting otions) */
            else if (IS_SACHDR_TIME(index)) {
                print_time(hd, index, modif, eyear, undef);
            }
            else if (IS_SACHDR_ENUM(index)) {
                print_enum(hd, index, modif, undef);
            }

            /* Process other float fields */
            else if (IS_SACHDR_FLOAT(index)) {
                fval = get_sac_hdr_float(&hd, index);
                if (fval != SAC_HEADER_FLOAT_UNDEFINED || !undef) {
                    printf("%g", fval);
                }
                else {
                    printf("%s", undef);
                }
            }

            /* Process other integer fields */
            else if (IS_SACHDR_INT(index)) {
                ival = get_sac_hdr_int(&hd, index);
                if (ival != SAC_HEADER_FLOAT_UNDEFINED || !undef) {
                    printf("%d", ival);
                }
                else {
                    printf("%s", undef);
                }
            }

            /* Process character fields */
            else {
                get_sac_hdr_string(&hd, index, rmblank, string);
                if (strncmp(string,SAC_HEADER_CHAR_UNDEFINED,6) || !undef) {
                    printf("%s", string);
                }
                else {
                    printf("%s", undef);
                }
            }
        } /* end of key-loop */


        /* Skip data part, if input from stdin (1 or 2 data sections) */
        if (!isfile) {
            ival = (!hd.leven || hd.iftype==IRLIM || hd.iftype==IAMPH) ? 2 : 1;
            for (j=0; j<ival*hd.npts; j++) {
                if (fread(&fval, sizeof(float), 1, stdin) != 1) {
                    printf("\n");
                    break;
                }
            }
        }

    } /* end of file-loop */


    /* Print newline */
    if (issac) printf("\n");


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



void print_enum(SACHEAD hd, int index, char modif, char *undef)
/************************************************************************
print_enum - print SAC enum header field value to stdout

Input:
hd      SAC header structure   (see SACHEAD typedef)
index   index of header field  (enum field assumed)
modif   output format modifier (see switch statement)
undef   string used for undefined header fields

Output: prints field value to stdout

Author: Nils Maercklin, 2008-10-05
*************************************************************************/
{
    int ival = get_sac_hdr_int(&hd, index);

    if (ival==SAC_HEADER_INT_UNDEFINED) {
        if (!undef) printf("%d", ival);
        else printf("%s", undef);
    }
    else {
        switch (modif) {
            case 'N':    /* Enum name */
            case 'n':
                printf("%s", sac_enum_name(ival, 0) );
                break;
            case 'D':    /* Enum description */
            case 'd':
                printf("%s", sac_enum_name(ival, 1) );
                break;
            case ' ':    /* Header field value */
            default:
                printf("%d", ival);
                break;
        }
    }
}



void print_time(SACHEAD hd, int index, char modif, int eyear, char *undef)
/************************************************************************
print_time - print SAC time header field value to stdout

Input:
hd      SAC header structure   (see SACHEAD typedef)
index   index of header field  (time), or INDEX_ZTIME
modif   output format modifier (see switch statement)
eyear   epoch year
undef   string used for undefined header fields

Output: prints field value to stdout

Author: Nils Maercklin, 2008-10-05
*************************************************************************/
{
    double etime=0.0;
    float fval=0.0;
    float sec=0.0;
    int month=0;
    int day=0;
    int year,hour,minute;
    int valid_jday=0;

    /* Get seconds of zerotime as float */
    sec = ((float)hd.nzsec) + 0.001*((float)hd.nzmsec);

    /* Convert Julian day to human date */
    if (check_zdate(hd)) {
        valid_jday = jday_to_hdate(hd.nzyear, hd.nzjday, &month, &day);
    }

    /* Get header value */
    if (index == INDEX_ZTIME) {
        fval = 0.0;
    }
    else {
        fval = get_sac_hdr_float(&hd, index);
    }


    /* Print time */
    if (fval == SAC_HEADER_FLOAT_UNDEFINED) {
        if (!undef) printf("%g", fval);
        else printf("%s", undef);   
    }
    else {
        switch (modif) {
            case 'E':    /* Epoch time in seconds */
            case 'e':
                if (check_zdate(hd) && check_ztime(hd)) {
                    etime = date_to_etime(eyear, hd.nzyear, hd.nzjday, \
                        hd.nzhour, hd.nzmin, sec) + ((double) fval);
                    printf("%.4lf", etime);
                }
                else {
                    if (!undef) printf("%g", fval);
                    else printf("%s", undef);   
                }
                break;
            case 'H':    /* Human date and time string */
            case 'h':
                if (valid_jday && check_ztime(hd)) {
                    etime = date_to_etime(eyear, hd.nzyear, hd.nzjday, \
                        hd.nzhour, hd.nzmin, sec) + ((double) fval);
                    etime_to_date(etime, eyear, &year, &day, &hour, \
                        &minute, &sec);
                    jday_to_hdate(year, day, &month, &day);
                    
                    printf("%04d-%02d-%02d %02d:%02d:%07.4f", \
                        year, month, day, hour, minute, sec);
                }
                else {
                    if (!undef) printf("%g", fval);
                    else printf("%s", undef);   
                }
                break;
            case 'N':    /* NORSAR time string */
            case 'n':
                if (check_zdate(hd) && check_ztime(hd)) {
                    etime = date_to_etime(eyear, hd.nzyear, hd.nzjday, \
                        hd.nzhour, hd.nzmin, sec) + ((double) fval);
                    etime_to_date(etime, eyear, &year, &day, &hour, \
                        &minute, &sec);
                    
                    printf("%04d-%03d:%02d.%02d.%06.3f", \
                        year, day, hour, minute, sec);
                }
                else {
                    if (!undef) printf("%g", fval);
                    else printf("%s", undef);   
                }
                break;
            case 'T':    /* Time of day (hour, minute, second) */
            case 't':
                if (check_ztime(hd)) {
                    etime = date_to_etime(eyear, eyear, 1, \
                        hd.nzhour, hd.nzmin, sec) + ((double) fval);
                    etime_to_date(etime, eyear, &year, &day, &hour, \
                        &minute, &sec);
                    
                    printf("%02d:%02d:%07.4f", hour, minute, sec);
                }
                else {
                    if (!undef) printf("%g", fval);
                    else printf("%s", undef);   
                }
                break;

            default:     /* Header value */
                printf("%g", fval);
                break;
        }
    }
}

/* END OF FILE */
