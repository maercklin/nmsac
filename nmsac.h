/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
NMSAC.H - Definitions of the SAC binary file header 
*************************************************************************

Author: Nils Maercklin,
  RISSC, University of Naples, Italy, September 2008

Version: 2010-01-30

Credits: These definitions are taken from "utils/sac.h" distributed  
  with SAC 101.2 (2008); "sac.h" was written by by Dennis O'Neill (1988), 
  and was later modified by Lorraine Hwang (1991), Xiaoming Ding (1993), 
  Lupei Zhu (1996), and Brian Savage (2008).

Modifications:
  2008-10-31 (NM): renamed fields xsize and ysiye to nxsize and nysize.
  2010-01-30 (NM): Moved global variables to "nmsacg.h".

Notes:
  This SAC header file is to be used with the routines in the NMSAC 
  software package. This file does not require the SAC distribution.
  The SAC header version is NVHDR = 6.

  Definitions should be compatible with routines in the "utils/" directory
  of the source distribution of SAC but may differ from those in the SacIO 
  library (SAC 101.2, 2008).

References:
  O'Neill, D. (1987). IRIS Interim Data Distribution Format (SAC ASCII), 
      Version 1.0 (12 November 1987). Incorporated Research Institutions 
      for Seismology (IRIS), 1616 North Fort Myer Drive, Suite 1440, 
      Arlington, Virginia 22209. 11 pages.
  Tull, J. (1987). SAC User's Manual, Version 10.2, October 7, 1987.
      Lawrence Livermore National Laboratory, L-205, Livermore, 
      California 94550.
  SAC - Seismic Analysis Code, User's manual, 
      http://www.iris.edu/manuals/sac/manual.html (last checked: 2008-02-27)


Comment flags describing each SAC header field (SACHEAD structure):
  Column 1:
      R = required by SAC
      (blank) = optional
  Column 2:
      A = settable from a priori knowledge
      D = available in data
      F = available in or derivable from SEED fixed data header
      T = available in SEED header tables
      (blank) = not directly available from SEED data, header, tables

*************************************************************************/

#ifndef _sachead_h
#define _sachead_h

/* True/false definitions */
#ifndef TRUE
#define FALSE   0
#define TRUE    1
#endif

/* Constants */
#define SAC_HEADER_FIELDS          133
#define SAC_HEADER_SIZE_NUMBERS    440
#define SAC_HEADER_SIZE            632
#define SAC_HEADER_TMARK_POSITION   10
#define SAC_HEADER_USERN_POSITION   40

#define SAC_HEADER_FLOAT_MIN         0
#define SAC_HEADER_FLOAT_MAX        69
#define SAC_HEADER_INT_MIN          70
#define SAC_HEADER_INT_MAX          84
#define SAC_HEADER_ENUM_MIN         85
#define SAC_HEADER_ENUM_MAX        104
#define SAC_HEADER_LOGICAL_MIN     105
#define SAC_HEADER_LOGICAL_MAX     109
#define SAC_HEADER_CHAR_MIN        110
#define SAC_HEADER_CHAR_MAX        133
#define SAC_HEADER_CHAR_DOUBLE     111
#define SAC_HEADER_CHAR_DOUBLE_END 112

#define SAC_HEADER_FLOAT_UNDEFINED (-12345.0)
#define SAC_HEADER_INT_UNDEFINED   (-12345)
#define SAC_HEADER_CHAR_UNDEFINED  ("-12345  ")
#define SAC_HEADER_UNDEFINED       ("UNDEF")

/* SAC header structure SACHEAD */
typedef struct sac_head
{
    float  delta;      /* RF time increment, sec    */
    float  depmin;     /*    minimum amplitude      */
    float  depmax;     /*    maximum amplitude      */
    float  scale;      /*    amplitude scale factor */
    float  odelta;     /*    observed time inc      */
    float  b;          /* RD initial time - wrt nz* */
    float  e;          /* RD end time               */
    float  o;          /*    event start            */
    float  a;          /*    1st arrival time       */
    float  fmt;        /*    internal use           */
    float  t0;         /*    user-defined time pick */
    float  t1;         /*    user-defined time pick */
    float  t2;         /*    user-defined time pick */
    float  t3;         /*    user-defined time pick */
    float  t4;         /*    user-defined time pick */
    float  t5;         /*    user-defined time pick */
    float  t6;         /*    user-defined time pick */
    float  t7;         /*    user-defined time pick */
    float  t8;         /*    user-defined time pick */
    float  t9;         /*    user-defined time pick */
    float  f;          /*    event end, sec > 0     */
    float  resp0;      /*    instrument respnse parm*/
    float  resp1;      /*    instrument respnse parm*/
    float  resp2;      /*    instrument respnse parm*/
    float  resp3;      /*    instrument respnse parm*/
    float  resp4;      /*    instrument respnse parm*/
    float  resp5;      /*    instrument respnse parm*/
    float  resp6;      /*    instrument respnse parm*/
    float  resp7;      /*    instrument respnse parm*/
    float  resp8;      /*    instrument respnse parm*/
    float  resp9;      /*    instrument respnse parm*/
    float  stla;       /*  T station latititude     */
    float  stlo;       /*  T station longitude      */
    float  stel;       /*  T station elevation, m   */
    float  stdp;       /*  T station depth, m       */
    float  evla;       /*    event latitude         */
    float  evlo;       /*    event longitude        */
    float  evel;       /*    event elevation        */
    float  evdp;       /*    event depth            */
    float  mag;        /*    magnitude value        */
    float  user0;      /*    available to user      */
    float  user1;      /*    available to user      */
    float  user2;      /*    available to user      */
    float  user3;      /*    available to user      */
    float  user4;      /*    available to user      */
    float  user5;      /*    available to user      */
    float  user6;      /*    available to user      */
    float  user7;      /*    available to user      */
    float  user8;      /*    available to user      */
    float  user9;      /*    available to user      */
    float  dist;       /*    stn-event distance, km */
    float  az;         /*    event-stn azimuth      */
    float  baz;        /*    stn-event azimuth      */
    float  gcarc;      /*    stn-event dist, degrees*/
    float  sb;         /*    saved b value          */
    float  sdelta;     /*    saved delta value      */
    float  depmen;     /*    mean value, amplitude  */
    float  cmpaz;      /*  T component azimuth      */
    float  cmpinc;     /*  T component inclination  */
    float  xminimum;   /*    XYZ X minimum value    */
    float  xmaximum;   /*    XYZ X maximum value    */
    float  yminimum;   /*    XYZ Y minimum value    */
    float  ymaximum;   /*    XYZ Y maximum value    */
    float  unused6;    /*    reserved for future use*/
    float  unused7;    /*    reserved for future use*/
    float  unused8;    /*    reserved for future use*/
    float  unused9;    /*    reserved for future use*/
    float  unused10;   /*    reserved for future use*/
    float  unused11;   /*    reserved for future use*/
    float  unused12;   /*    reserved for future use*/
    int    nzyear;     /*  F zero time of file, yr  */
    int    nzjday;     /*  F zero time of file, day */
    int    nzhour;     /*  F zero time of file, hr  */
    int    nzmin;      /*  F zero time of file, min */
    int    nzsec;      /*  F zero time of file, sec */
    int    nzmsec;     /*  F zero time of file, msec*/
    int    nvhdr;      /*  R header version number  */
    int    norid;      /*    Origin ID              */
    int    nevid;      /*    Event ID               */
    int    npts;       /* RF number of samples      */
    int    nsnpts;     /*    saved npts             */
    int    nwfid;      /*    Waveform ID            */
    int    nxsize;     /*    XYZ X size             */
    int    nysize;     /*    XYZ Y size             */
    int    unused15;   /*    reserved for future use*/
    int    iftype;     /* RA type of file           */
    int    idep;       /*    type of amplitude      */
    int    iztype;     /*    zero time equivalence  */
    int    unused16;   /*    reserved for future use*/
    int    iinst;      /*    recording instrument   */
    int    istreg;     /*    stn geographic region  */
    int    ievreg;     /*    event geographic region*/
    int    ievtyp;     /*    event type             */
    int    iqual;      /*    quality of data        */
    int    isynth;     /*    synthetic data flag    */
    int    imagtyp;    /*    magnitude type         */
    int    imagsrc;    /*    magnitude source       */
    int    unused19;   /*    reserved for future use*/
    int    unused20;   /*    reserved for future use*/
    int    unused21;   /*    reserved for future use*/
    int    unused22;   /*    reserved for future use*/
    int    unused23;   /*    reserved for future use*/
    int    unused24;   /*    reserved for future use*/
    int    unused25;   /*    reserved for future use*/
    int    unused26;   /*    reserved for future use*/
    int    leven;      /* RA data-evenly-spaced flag*/
    int    lpspol;     /*    station polarity flag  */
    int    lovrok;     /*    overwrite permission   */
    int    lcalda;     /*    calc distance, azimuth */
    int    unused27;   /*    reserved for future use*/
    char   kstnm[8];   /*  F station name           */
    char   kevnm[16];  /*    event name             */
    char   khole[8];   /*    man-made event name    */
    char   ko[8];      /*    event origin time id   */
    char   ka[8];      /*    1st arrival time ident */
    char   kt0[8];     /*    time pick 0 ident      */
    char   kt1[8];     /*    time pick 1 ident      */
    char   kt2[8];     /*    time pick 2 ident      */
    char   kt3[8];     /*    time pick 3 ident      */
    char   kt4[8];     /*    time pick 4 ident      */
    char   kt5[8];     /*    time pick 5 ident      */
    char   kt6[8];     /*    time pick 6 ident      */
    char   kt7[8];     /*    time pick 7 ident      */
    char   kt8[8];     /*    time pick 8 ident      */
    char   kt9[8];     /*    time pick 9 ident      */
    char   kf[8];      /*    end of event ident     */
    char   kuser0[8];  /*    available to user      */
    char   kuser1[8];  /*    available to user      */
    char   kuser2[8];  /*    available to user      */
    char   kcmpnm[8];  /*  F component name         */
    char   knetwk[8];  /*    network name           */
    char   kdatrd[8];  /*    date data read         */
    char   kinst[8];   /*    instrument name        */
} SACHEAD;


/* Enumerated SAC header values */
enum SAC_HEADER_ENUMS {
    /* enumerated header values */
    IREAL    = 0,   /* Undocumented                */
    ITIME    = 1,   /* Time series file            */
    IRLIM    = 2,   /* Spectral file-real/imag     */
    IAMPH    = 3,   /* Spectral file-ampl/phase    */
    IXY      = 4,   /* General x vs y file         */
    IUNKN    = 5,   /* Unknown                     */
    IDISP    = 6,   /* Displacement (NM)           */
    IVEL     = 7,   /* Velocity (NM/SEC)           */
    IACC     = 8,   /* Acceleration (CM/SEC/SEC)   */
    IB       = 9,   /* Begin time                  */
    IDAY     = 10,  /* GMT day                     */
    IO       = 11,  /* Event origin time           */
    IA       = 12,  /* First arrival time          */
    IT0      = 13,  /* User defined time pick 0    */
    IT1      = 14,  /* User defined time pick 1    */
    IT2      = 15,  /* User defined time pick 2    */
    IT3      = 16,  /* User defined time pick 3    */
    IT4      = 17,  /* User defined time pick 4    */
    IT5      = 18,  /* User defined time pick 5    */
    IT6      = 19,  /* User defined time pick 6    */
    IT7      = 20,  /* User defined time pick 7    */
    IT8      = 21,  /* User defined time pick 8    */
    IT9      = 22,  /* User defined time pick 9    */
    IRADNV   = 23,  /* Radial (NTS)                */
    ITANNV   = 24,  /* Tangential (NTS)            */
    IRADEV   = 25,  /* Radial (EVENT)              */
    ITANEV   = 26,  /* Tangential (EVENT)          */
    INORTH   = 27,  /* North positive              */
    IEAST    = 28,  /* East positive               */
    IHORZA   = 29,  /* Horizontal (ARB)            */
    IDOWN    = 30,  /* Down positive               */
    IUP      = 31,  /* Up positive                 */
    ILLLBB   = 32,  /* LLL broadband               */
    IWWSN1   = 33,  /* WWSN 15-100                 */
    IWWSN2   = 34,  /* WWSN 30-100                 */
    IHGLP    = 35,  /* High-gain long-period       */
    ISRO     = 36,  /* SRO                         */
    INUCL    = 37,  /* Nuclear event               */
    IPREN    = 38,  /* Nuclear pre-shot event      */
    IPOSTN   = 39,  /* Nuclear post-shot event     */
    IQUAKE   = 40,  /* Earthquake                  */
    IPREQ    = 41,  /* Foreshock                   */
    IPOSTQ   = 42,  /* Aftershock                  */
    ICHEM    = 43,  /* Chemical explosion          */
    IOTHER   = 44,  /* Other                       */
    IGOOD    = 45,  /* Good                        */
    IGLCH    = 46,  /* Gliches                     */
    IDROP    = 47,  /* Dropouts                    */
    ILOWSN   = 48,  /* Low signal to noise ratio   */
    IRLDTA   = 49,  /* Real data                   */
    IVOLTS   = 50,  /* Velocity (volts)            */
    IXYZ     = 51,  /* General XYZ (3-D) file      */
    /* These 18 added to describe magnitude type and source maf 970205 */
    IMB      = 52,  /* Bodywave Magnitude          */
    IMS      = 53,  /* Surface Magnitude           */
    IML      = 54,  /* Local Magnitude             */
    IMW      = 55,  /* Moment Magnitude            */
    IMD      = 56,  /* Duration Magnitude          */
    IMX      = 57,  /* User Defined Magnitude      */
    INEIC    = 58,  /* INEIC                       */
    IPDEQ    = 59,  /* IPDEQ                       */
    IPDEW    = 60,  /* IPDEW                       */
    IPDE     = 61,  /* IPDE                        */
    IISC     = 62,  /* IISC                        */
    IREB     = 63,  /* IREB                        */
    IUSGS    = 64,  /* IUSGS                       */
    IBRK     = 65,  /* IBRK                        */
    ICALTECH = 66,  /* ICALTECH                    */
    ILLNL    = 67,  /* ILLNL                       */
    IEVLOC   = 68,  /* IEVLOC                      */
    IJSOP    = 69,  /* IJSOP                       */
    IUSER    = 70,  /* IUSER                       */
    IUNKNOWN = 71,  /* IUNKNOWN                    */
    /* These  17 added for ievtyp. maf 970325 */
    IQB      = 72,  /* Quarry or mine blast confirmed by quarry */
    IQB1     = 73,  /* Quarry or mine blast with designed shot information-ripple fired */
    IQB2     = 74,  /* Quarry or mine blast with observed shot information-ripple fired */
    IQBX     = 75,  /* Quarry or mine blast - single shot */
    IQMT     = 76,  /* Quarry or mining-induced events: tremors and rockbursts */
    IEQ      = 77,  /* Earthquake                  */
    IEQ1     = 78,  /* Earthquakes in a swarm or aftershock sequence */
    IEQ2     = 79,  /* Felt earthquake             */
    IME      = 80,  /* Marine explosion            */
    IEX      = 81,  /* Other explosion             */
    INU      = 82,  /* Nuclear explosion           */
    INC      = 83,  /* Nuclear cavity collapse     */
    IO_      = 84,  /* Other source of known origin */
    IL       = 85,  /* Local event of unknown origin */
    IR       = 86,  /* Regional event of unknown origin */
    IT       = 87,  /* Teleseismic event of unknown origin */
    IU       = 88,  /* Undetermined or conflicting information  */
    /* These 9 added for ievtype to keep up with database. maf 000530 */
    IEQ3     = 89,  /* Damaging Earthquake         */
    IEQ0     = 90,  /* Probable earthquake         */
    IEX0     = 91,  /* Probable explosion          */
    IQC      = 92,  /* Mine collapse               */
    IQB0     = 93,  /* Probable Mine Blast         */
    IGEY     = 94,  /* Geyser                      */
    ILIT     = 95,  /* Light                       */
    IMET     = 96,  /* Meteroic event              */
    IODOR    = 97   /* Odors                       */
};


enum SAC_HEADER_TYPES {
  SAC_HEADER_UNDEFINED_TYPE = 1, 
  SAC_HEADER_FLOAT_TYPE,   SAC_HEADER_INT_TYPE,   SAC_HEADER_ENUM_TYPE,
  SAC_HEADER_LOGICAL_TYPE, SAC_HEADER_CHAR8_TYPE, SAC_HEADER_CHAR16_TYPE
};


enum SAC_HEADER_ORDER {
  SAC_HEADER_DELTA = 0, 
  SAC_HEADER_DEPMIN,   SAC_HEADER_DEPMAX,    SAC_HEADER_SCALE,    SAC_HEADER_ODELTA,   SAC_HEADER_B,
  SAC_HEADER_E,        SAC_HEADER_O,         SAC_HEADER_A,        SAC_HEADER_FMT,      SAC_HEADER_T0,       
  SAC_HEADER_T1,       SAC_HEADER_T2,        SAC_HEADER_T3,       SAC_HEADER_T4,       SAC_HEADER_T5,       
  SAC_HEADER_T6,       SAC_HEADER_T7,        SAC_HEADER_T8,       SAC_HEADER_T9,       SAC_HEADER_F,        
  SAC_HEADER_RESP0,    SAC_HEADER_RESP1,     SAC_HEADER_RESP2,    SAC_HEADER_RESP3,    SAC_HEADER_RESP4,    
  SAC_HEADER_RESP5,    SAC_HEADER_RESP6,     SAC_HEADER_RESP7,    SAC_HEADER_RESP8,    SAC_HEADER_RESP9,    
  SAC_HEADER_STLA,     SAC_HEADER_STLO,      SAC_HEADER_STEL,     SAC_HEADER_STDP,     SAC_HEADER_EVLA,     
  SAC_HEADER_EVLO,     SAC_HEADER_EVEL,      SAC_HEADER_EVDP,     SAC_HEADER_MAG,      SAC_HEADER_USER0,    
  SAC_HEADER_USER1,    SAC_HEADER_USER2,     SAC_HEADER_USER3,    SAC_HEADER_USER4,    SAC_HEADER_USER5,    
  SAC_HEADER_USER6,    SAC_HEADER_USER7,     SAC_HEADER_USER8,    SAC_HEADER_USER9,    SAC_HEADER_DIST,     
  SAC_HEADER_AZ,       SAC_HEADER_BAZ,       SAC_HEADER_GCARC,    SAC_HEADER_SB,       SAC_HEADER_SDELTA,  
  SAC_HEADER_DEPMEN,   SAC_HEADER_CMPAZ,     SAC_HEADER_CMPINC,   SAC_HEADER_XMINIMUM, SAC_HEADER_XMAXIMUM, 
  SAC_HEADER_YMINIMUM, SAC_HEADER_YMAXIMUM,  SAC_HEADER_UNUSED6,  SAC_HEADER_UNUSED7,  SAC_HEADER_UNUSED8,  
  SAC_HEADER_UNUSED9,  SAC_HEADER_UNUSED10,  SAC_HEADER_UNUSED11, SAC_HEADER_UNUSED12, SAC_HEADER_NZYEAR,
  SAC_HEADER_NZJDAY,   SAC_HEADER_NZHOUR,    SAC_HEADER_NZMIN,    SAC_HEADER_NZSEC,    SAC_HEADER_NSMSEC,
  SAC_HEADER_NVHDR,    SAC_HEADER_NORID,     SAC_HEADER_NEVID,    SAC_HEADER_NPTS,     SAC_HEADER_NSNPTS,
  SAC_HEADER_NWFID,    SAC_HEADER_NXSIZE,    SAC_HEADER_NYSIZE,   SAC_HEADER_UNUSED15, SAC_HEADER_IFTYPE,
  SAC_HEADER_IDEP,     SAC_HEADER_IZTYPE,    SAC_HEADER_UNUSED16, SAC_HEADER_IINST,    SAC_HEADER_ISTREG, 
  SAC_HEADER_IEVREG,   SAC_HEADER_IEVTYP,    SAC_HEADER_IQUAL,    SAC_HEADER_ISYNTH,   SAC_HEADER_IMAGTYPE, 
  SAC_HEADER_IMAGSRC,  SAC_HEADER_UNUSED19,  SAC_HEADER_UNUSED20, SAC_HEADER_UNUSED21, SAC_HEADER_UNUSED22, 
  SAC_HEADER_UNUSED23, SAC_HEADER_UNUSED24,  SAC_HEADER_UNUSED25, SAC_HEADER_UNUSED26, SAC_HEADER_LEVEN,
  SAC_HEADER_LPSPOL,   SAC_HEADER_LOVROK,    SAC_HEADER_LCALDA,   SAC_HEADER_UNUSED27, SAC_HEADER_KSTNM,
  SAC_HEADER_KEVNM,    SAC_HEADER_KEVNM_END, SAC_HEADER_KHOLE,    SAC_HEADER_KO,       SAC_HEADER_KA,
  SAC_HEADER_KT0,      SAC_HEADER_KT1,       SAC_HEADER_KT2,      SAC_HEADER_KT3,      SAC_HEADER_KT4,
  SAC_HEADER_KT5,      SAC_HEADER_KT6,       SAC_HEADER_KT7,      SAC_HEADER_KT8,      SAC_HEADER_KT9,
  SAC_HEADER_KF,       SAC_HEADER_KUSER0,    SAC_HEADER_KUSER1,   SAC_HEADER_KUSER2,   SAC_HEADER_KCMPNM,
  SAC_HEADER_KNETWK,   SAC_HEADER_KDATRD,    SAC_HEADER_KINST
};

#endif

/* END OF FILE */
