/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
NMETIME.C - Functions to calculate epoch time from human date and time, 
            and vice-versa
*************************************************************************

Author: Nils Maercklin,
  RISSC, University of Naples, Italy, October 2008

Version: 2008-10-05

*************************************************************************/

/* Include files */
#include <stdlib.h>
#include <math.h>


/* Function prototypes */
int ndays_yearmonth(int year, int month);
int hdate_to_jday(int year, int month, int day);
int jday_to_hdate(int year, int jday, int *month, int *day);
double date_to_etime(int eyear, int year, int jday, int hour, int minute, \
    float sec);
void etime_to_date(double etime, int eyear, int *year, int *jday, \
    int *hour, int *minute, float *sec);


/************************************************************************/
/* Functions for epoch time conversion                                  */
/************************************************************************/

int ndays_yearmonth(int year, int month)
/************************************************************************
ndays_yearmonth - check for leap year and return days in year or month

Input:
year    year
month   month (1-12; or 0 to get number of days of specified year)

Output: returns the number of days (0 in case of invalid month)

Author: Nils Maercklin, 2004-07-19
*************************************************************************/
{
    const int ndays[13] = {365,31,28,31,30,31,30,31,31,30,31,30,31};
    int leap = (((year%4 == 0 && year%100 != 0) || (year%400 == 0)) ? 1 : 0);

    if (month<0  || month>12) return 0;
    if (month==0 || month==2) return ndays[month] + leap;

    return ndays[month];
}



int hdate_to_jday(int year, int month, int day)
/************************************************************************
hdate_to_jday - get day of the year (Julian) from month and day

Input:
year    year
month   month (1-12)
day     day of month (positive number)

Output: returns Julian day (1-366, 0 in case of invalid month)

Author: Nils Maercklin, 2004-07-19
*************************************************************************/
{
    int i, jday;

    if (month<1 || month>12) return 0;

    jday = day;
    for (i=1; i<month; i++) jday += ndays_yearmonth(year, i);

    return jday;
}



int jday_to_hdate(int year, int jday, int *month, int *day)
/************************************************************************
jday_to_hdate - get day and month from day of the year (Julian)

Input:
year    year
jday    Julian day (1-365 or 1-366)

Output: returns 1 (success) or 0 (invalid jday)
month   month (1-12)
day     day of month (1-31)

Author: Nils Maercklin, 2004-11-14
*************************************************************************/
{
    int i;

    if (jday<1 || jday>ndays_yearmonth(year,0)) {
        *month = 0;
        *day   = 0;
        return 0;
    }

    for (i=1; i<=12; i++) {
        *month = i;
        if (jday <= ndays_yearmonth(year,i)) break;
        jday -= ndays_yearmonth(year,i);
    }

    *day = jday;

    return 1;
}



double date_to_etime(int eyear, int year, int jday, int hour, int minute, \
    float sec)
/************************************************************************
date_to_etime - calculate epoch time from date and time of day

Input:
eyear   epoch year, e.g. 1970
year    year
jday    Julian day (1-365 or 1-366)
hour    hour of day (24 hour clock)
minute  minutes of hour
sec     seconds

Output: returns epoch time in seconds

Author: Nils Maercklin, 2008-10-03
*************************************************************************/
{
    int i;
    int edays=0;       /* number of days since epoch start */
    double etime=0.0;  /* epoch time in seconds */

    /* Epoch days */
    if (year>=eyear) {
        for (i=eyear; i<year; i++) edays += ndays_yearmonth(i,0);
    }
    else {
        for (i=year; i<eyear; i++) edays -= ndays_yearmonth(i,0);
    }
    edays += jday-1;

    /* Convert epoch days and time of day to epoch time */
    etime  = 86400.0 * ((double)edays);
    etime += 3600.0*((double)hour) + 60.0*((double)minute) + ((double)sec);

    return etime;
}



void etime_to_date(double etime, int eyear, int *year, int *jday, \
    int *hour, int *minute, float *sec)
/************************************************************************
etime_to_date - convert epoch time to date and time of day

Input:
etime   epoch time in seconds
eyear   epoch year, e.g. 1970

Output:
year    year
jday    Julian day (1-365 or 1-366)
hour    hour of day (24 hour clock)
minute  minutes of hour
sec     seconds

Author: Nils Maercklin, 2008-10-03
*************************************************************************/
{
    int edays=0;       /* number of days since epoch start */
    int ndays=0;       /* number of days of a given year */
    double daysec=0.0; /* seconds of the last day */

    /* Epoch days and remaining seconds of last day */
    edays  = (int) floor(etime/86400.0);
    daysec = fmod(etime, 86400.0);

    /* Get year */
    if (edays>=0) {
        while ((ndays=ndays_yearmonth(eyear,0)) <= edays) {
            edays -= ndays;
            eyear++;
        }
    }
    else {
        if (daysec<0.0) daysec += 86400.0;
        do {
            eyear--;
            ndays  = ndays_yearmonth(eyear,0);
            edays += ndays;
        } while (edays<0);
    }

    /* Year and Julian day */
    *year   = eyear;
    *jday   = edays + 1;

    /* Time of day */
    *hour   = (int) floor(daysec/3600.0);
    daysec -= 3600.0*floor(daysec/3600.0);
    *minute = (int) floor(daysec/60.0);
    daysec -= 60.0*floor(daysec/60.0);
    *sec    = (float) daysec;
}

/* END OF FILE */
