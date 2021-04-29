/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
NMSACNAM.H - Global variables with names of SAC header fields and values
*************************************************************************

Author: Nils Maercklin,
  RISSC, University of Naples, Italy, January 2010

Version: 2010-01-30

Credits: These variables are taken from "utils/sac.h" distributed  
  with SAC 101.2 (2008); "sac.h" was written by by Dennis O'Neill (1988), 
  and was later modified by Lorraine Hwang (1991), Xiaoming Ding (1993), 
  Lupei Zhu (1996), and Brian Savage (2008).

Notes:
  This file is to be used only with the NMSAC library functions defined
  in the file "nmsaclib.c". The SAC header version is NVHDR = 6.

  These variables should be compatible with routines in the "utils/" 
  directory of the source distribution of SAC (SAC 101.2, 2008).
  

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

*************************************************************************/


char *SacHeaderNameNull = "";

/* SAC header field names */
char *SacHeaderName[] = {
    /* floats */
    "delta",       /* RF time increment, sec    */
    "depmin",      /*    minimum amplitude      */
    "depmax",      /*    maximum amplitude      */
    "scale",       /*    amplitude scale factor */
    "odelta",      /*    observed time inc      */
    "b",           /* RD initial time - wrt nz* */
    "e",           /* RD end time               */
    "o",           /*    event start            */
    "a",           /*    1st arrival time       */
    "Fmt",         /*    internal use           */
    "t0",          /*    user-defined time pick */
    "t1",          /*    user-defined time pick */
    "t2",          /*    user-defined time pick */
    "t3",          /*    user-defined time pick */
    "t4",          /*    user-defined time pick */
    "t5",          /*    user-defined time pick */
    "t6",          /*    user-defined time pick */
    "t7",          /*    user-defined time pick */
    "t8",          /*    user-defined time pick */
    "t9",          /*    user-defined time pick */
    "F",           /*    event end, sec > 0     */
    "resp0",       /*    instrument respnse parm*/
    "resp1",       /*    instrument respnse parm*/
    "resp2",       /*    instrument respnse parm*/
    "resp3",       /*    instrument respnse parm*/
    "resp4",       /*    instrument respnse parm*/
    "resp5",       /*    instrument respnse parm*/
    "resp6",       /*    instrument respnse parm*/
    "resp7",       /*    instrument respnse parm*/
    "resp8",       /*    instrument respnse parm*/
    "resp9",       /*    instrument respnse parm*/
    "stla",        /*  T station latititude     */
    "stlo",        /*  T station longitude      */
    "stel",        /*  T station elevation, m   */
    "stdp",        /*  T station depth, m       */
    "evla",        /*    event latitude         */
    "evlo",        /*    event longitude        */
    "evel",        /*    event elevation        */
    "evdp",        /*    event depth            */
    "mag",         /*    reserved for future use*/
    "user0",       /*    available to user      */
    "user1",       /*    available to user      */
    "user2",       /*    available to user      */
    "user3",       /*    available to user      */
    "user4",       /*    available to user      */
    "user5",       /*    available to user      */
    "user6",       /*    available to user      */
    "user7",       /*    available to user      */
    "user8",       /*    available to user      */
    "user9",       /*    available to user      */
    "dist",        /*    stn-event distance, km */
    "az",          /*    event-stn azimuth      */
    "baz",         /*    stn-event azimuth      */
    "gcarc",       /*    stn-event dist, degrees*/
    "sb",          /*    internal use           */
    "sdelta",      /*    internal use           */
    "depmen",      /*    mean value, amplitude  */
    "cmpaz",       /*  T component azimuth      */
    "cmpinc",      /*  T component inclination  */
    "xminimum",    /*    reserved for future use*/
    "xmaximum",    /*    reserved for future use*/
    "yminimum",    /*    reserved for future use*/
    "ymaximum",    /*    reserved for future use*/
    "unused6",     /*    reserved for future use*/
    "unused7",     /*    reserved for future use*/
    "unused8",     /*    reserved for future use*/
    "unused9",     /*    reserved for future use*/
    "unused10",    /*    reserved for future use*/
    "unused11",    /*    reserved for future use*/
    "unused12",    /*    reserved for future use*/  
    /* ints */
    "nzyear",      /*  F zero time of file, yr  */
    "nzjday",      /*  F zero time of file, day */
    "nzhour",      /*  F zero time of file, hr  */
    "nzmin",       /*  F zero time of file, min */
    "nzsec",       /*  F zero time of file, sec */
    "nzmsec",      /*  F zero time of file, msec*/
    "nvhdr",       /*  R header version number  */
    "norid",       /*    internal use           */
    "nevid",       /*    internal use           */
    "npts",        /* RF number of samples      */
    "nsnpts",      /*    internal use           */
    "nwfid",       /*    internal use           */
    "nxsize",      /*    reserved for future use  (NM: xsize => nxsize) */
    "nysize",      /*    reserved for future use  (NM: ysize => nysize) */
    "unused15",    /*    reserved for future use*/
    "iftype",      /* RA type of file           */
    "idep",        /*    type of amplitude      */
    "iztype",      /*    zero time equivalence  */
    "unused16",    /*    reserved for future use*/
    "iinst",       /*    recording instrument   */
    "istreg",      /*    stn geographic region  */
    "ievreg",      /*    event geographic region*/
    "ievtyp",      /*    event type             */
    "iqual",       /*    quality of data        */
    "isynth",      /*    synthetic data flag    */
    "imagtyp",     /*    reserved for future use*/
    "imagsrc",     /*    reserved for future use*/
    "unused19",    /*    reserved for future use*/
    "unused20",    /*    reserved for future use*/
    "unused21",    /*    reserved for future use*/
    "unused22",    /*    reserved for future use*/
    "unused23",    /*    reserved for future use*/
    "unused24",    /*    reserved for future use*/
    "unused25",    /*    reserved for future use*/
    "unused26",    /*    reserved for future use*/
    "leven",       /* RA data-evenly-spaced flag*/
    "lpspol",      /*    station polarity flag  */
    "lovrok",      /*    overwrite permission   */
    "lcalda",      /*    calc distance, azimuth */
    "unused27",    /*    reserved for future use*/
    /* chars */
    "kstnm",       /*  F station name           */
    "kevnm",       /*    event name             */
    "kevnm empty", /*                           */
    "khole",       /*    man-made event name    */
    "ko",          /*    event origin time id   */
    "ka",          /*    1st arrival time ident */
    "kt0",         /*    time pick 0 ident      */
    "kt1",         /*    time pick 1 ident      */
    "kt2",         /*    time pick 2 ident      */
    "kt3",         /*    time pick 3 ident      */
    "kt4",         /*    time pick 4 ident      */
    "kt5",         /*    time pick 5 ident      */
    "kt6",         /*    time pick 6 ident      */
    "kt7",         /*    time pick 7 ident      */
    "kt8",         /*    time pick 8 ident      */
    "kt9",         /*    time pick 9 ident      */
    "kf",          /*    end of event ident     */
    "kuser0",      /*    available to user      */
    "kuser1",      /*    available to user      */
    "kuser2",      /*    available to user      */
    "kcmpnm",      /*  F component name         */
    "knetwk",      /*    network name           */
    "kdatrd",      /*    date data read         */
    "kinst"        /*    instrument name        */
};


/* SAC header enum names */
char *SacHeaderEnums[] = {
    "IREAL",       /* 0    To be consistent with defines above */
    /* iftype */
    "ITIME",       /* 1    Time series file            */
    "IRLIM",       /* 2    Spectral file-real/imag     */
    "IAMPH",       /* 3    Spectral file-ampl/phase    */
    "IXY",         /* 4    General x vs y file         */
    "IUNKN",       /* 5    Unknown                     */
    /* idep */
    "IDISP",       /* 6    Displacement (NM)           */
    "IVEL",        /* 7    Velocity (NM/SEC)           */
    "IACC",        /* 8    Acceleration (CM/SEC/SEC)   */
    /* iztype */
    "IB",          /* 9    Begin time                  */
    "IDAY",        /* 10   GMT day                     */
    "IO",          /* 11   Event origin time           */
    "IA",          /* 12   First arrival time          */
    "IT0",         /* 13   User defined time pick 0    */
    "IT1",         /* 14   User defined time pick 1    */
    "IT2",         /* 15   User defined time pick 2    */
    "IT3",         /* 16   User defined time pick 3    */
    "IT4",         /* 17   User defined time pick 4    */
    "IT5",         /* 18   User defined time pick 5    */
    "IT6",         /* 19   User defined time pick 6    */
    "IT7",         /* 20   User defined time pick 7    */
    "IT8",         /* 21   User defined time pick 8    */
    "IT9",         /* 22   User defined time pick 9    */
    /* iinst */
    "IRADNV",      /* 23   Radial (NTS)                */
    "ITANNV",      /* 24   Tangential (NTS)            */
    "IRADEV",      /* 25   Radial (EVENT)              */
    "ITANEV",      /* 26   Tangential (EVENT)          */
    "INORTH",      /* 27   North positive              */
    "IEAST",       /* 28   East positive               */
    "IHORZA",      /* 29   Horizontal (ARB)            */
    "IDOWN",       /* 30   Down positive               */
    "IUP",         /* 31   Up positive                 */
    "ILLLBB",      /* 32   LLL broadband               */
    "IWWSN1",      /* 33   WWSN 15-100                 */
    "IWWSN2",      /* 34   WWSN 30-100                 */
    "IHGLP",       /* 35   High-gain long-period       */
    "ISRO",        /* 36   SRO                         */
    /* ievtyp */
    "INUCL",       /* 37   Nuclear event               */
    "IPREN",       /* 38   Nuclear pre-shot event      */
    "IPOSTN",      /* 39   Nuclear post-shot event     */
    "IQUAKE",      /* 40   Earthquake                  */
    "IPREQ",       /* 41   Foreshock                   */
    "IPOSTQ",      /* 42   Aftershock                  */
    "ICHEM",       /* 43   Chemical explosion          */
    "IOTHER",      /* 44   Other                       */
    /* iqual */
    "IGOOD",       /* 45   Good                        */
    "IGLCH",       /* 46   Gliches                     */
    "IDROP",       /* 47   Dropouts                    */
    "ILOWSN",      /* 48   Low signal to noise ratio   */
    /* isynth */
    "IRLDTA",      /* 49   Real data                   */
    "IVOLTS",      /* 50   Velocity (volts)            */
    "IXYZ",        /* 51   General XYZ (3-D) file      */
    /* These 18 added to describe magnitude type and source maf 970205 */
    "IMB",         /* 52   Bodywave Magnitude */
    "IMS",         /* 53   Surface Magnitude */
    "IML",         /* 54   Local Magnitude */
    "IMW",         /* 55   Moment Magnitude */
    "IMD",         /* 56   Duration Magnitude */
    "IMX",         /* 57   User Defined Magnitude */
    "INEIC",       /* 58   INEIC */
    "IPDEQ",       /* 59   IPDEQ */
    "IPDEW",       /* 60   IPDEW */
    "IPDE",        /* 61   IPDE */
    "IISC",        /* 62   IISC */
    "IREB",        /* 63   IREB */
    "IUSGS",       /* 64   IUSGS */
    "IBRK",        /* 65   IBRK */
    "ICALTECH",    /* 66   ICALTECH */
    "ILLNL",       /* 67   ILLNL */
    "IEVLOC",      /* 68   IEVLOC */
    "IJSOP",       /* 69   IJSOP */
    "IUSER",       /* 70   IUSER */
    "IUNKNOWN",    /* 71   IUNKNOWN */
    /* These 17 added for ievtyp. maf 970325 */
    "IQB",         /* 72   Quarry or mine blast confirmed by quarry */
    "IQB1",        /* 73   Quarry or mine blast with designed shot information-ripple fired*/
    "IQB2",        /* 74   Quarry or mine blast with observed shot information-ripple fired*/
    "IQBX",        /* 75   Quarry or mine blast - single shot */
    "IQMT",        /* 76   Quarry or mining-induced events: tremors and rockbursts */
    "IEQ",         /* 77   Earthquake */
    "IEQ1",        /* 78   Earthquakes in a swarm or aftershock sequence */
    "IEQ2",        /* 79   Felt earthquake */
    "IME",         /* 80   Marine explosion */
    "IEX",         /* 81   Other explosion */
    "INU",         /* 82   Nuclear explosion */
    "INC",         /* 83   Nuclear cavity collapse */
    "IO_",         /* 84   Other source of known origin */
    "IL",          /* 85   Local event of unknown origin */
    "IR",          /* 86   Regional event of unknown origin */
    "IT",          /* 87   Teleseismic event of unknown origin */
    "IU",          /* 88   Undetermined or conflicting information  */
    /* These 9 added for ievtype to keep up with database. maf 000530 */
    "IEQ3",        /* 89   Damaging Earthquake */
    "IEQ0",        /* 90   Probable earthquake */
    "IEX0",        /* 91   Probable explosion */
    "IQC",         /* 92   Mine collapse */
    "IQB0",        /* 93   Probable Mine Blast */
    "IGEY",        /* 94   Geyser */
    "ILIT",        /* 95   Light */
    "IMET",        /* 96   Meteroic event */
    "IODOR"        /* 97   Odors */
};


/* Enum size (maximum enum value) */
const int SacHeaderEnumsLength = sizeof(SacHeaderEnums) / sizeof(char *);


/* SAC header enum descriptions */
char *SacHeaderEnumsDescription[] = {
    "Undocumented", 
    /* iftype */
    "Time Series File", "Spectral File-Real/Imag", "Spectral File-Ampl/Phase",
    "General X vs Y file", "Unknown", 
    /* idep */
    "Displacement (nm)", "Velocity (nm/sec)", "Acceleration (cm/sec/sec)",
    /* iztype */
    "Begin Time", "GMT Day", "Event Origin Time", "First Arrival Time",
    "User Defined Time Pick 0", "User Defined Time Pick 1",
    "User Defined Time Pick 2", "User Defined Time Pick 3",
    "User Defined Time Pick 4", "User Defined Time Pick 5",
    "User Defined Time Pick 6", "User Defined Time Pick 7",
    "User Defined Time Pick 8", "User Defined Time Pick 9",
    /* iinst */
    "Radial (NTS)", "Tangential (NTS)", "Radial (Event)", "Tangential (Event)",
    "North Positive", "East Positive", "Horizontal (ARB)", "Down Positive",
    "Up Positive", "LLL Broadband", "WWSN 15-100", "WWSN 30-100",
    "High Gain Long Period", "SRO",
    /* ievtyp */
    "Nuclear Event", "Nuclear Pre-Shot Event", "Nuclear Post-Shot Event",
    "Earthquake", "Foreshock", "Aftershock", "Chemical Explosion",
    "Other",
     /* iqual */
    "Good", "Glitches", "Dropouts", "Low Signal to Noise Ratio",
     /* isynth */
    "Real Data", "Velocity (Volts)", "General XYZ (3-D) file",
     /* These 18 added to describe magnitude type and source maf 970205 */
    "Body Wave Magnitude (mb)", "Surface Wave Magnitude (Ms)",
    "Local Magnitude (ML)", "Moment Magnitude (Mw)",
    "Duration Magnitude (Md)", "User Defined Magnitude",
    "NEIC", "PDEQ", "PDEW", "PDE", "ISC", "REB", "USGS", "Berkeley",
    "Caltech", "LLNL", "EVLOC", "JSOP", "User", "Unknown",
    /* These 17 added for ievtyp. maf 970325 */
    "Quarry/Mine Blast, Confirmed by Quarry", 
    "Quarry/Mine Blast with Shot Information, Ripple Fired", 
    "Quarry/Mine Blast with Observed Shot Information, Ripple Fired",
    "Quarry/Mine Blast, Single Shot",
    "Quarry or Mining Induced Events, Tremors and Rockbursts",
    "Earthquake", "Earthquake, Swarm or Aftershock Sequence",
    "Earthquake, Felt", "Marine Explosion",
    "Other Explosion", "Nuclear Explosion",
    "Nuclear Cavity Collapse", "Other Source, Unknown Origin",
    "Local Event, Unknown Origin", "Regional Event, Unknown Origin",
    "Teleseismic Event, Unknown Origin", 
    "Undetermined or Conflicting Information",
    /* These 9 added for ievtype to keep up with database. maf 000530 */
    "Damaging Earthquake", "Probable Earthquake", "Probable Explosion",
    "Mine Collapse", "Probable Mine Blast", "Geyser",    
    "Light", "Meteroic Event", "Odors" 
};

/* END OF FILE */
