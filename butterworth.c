/*****************************************************************************
BUTTERWORTH - Functions to design and apply Butterworth filters
******************************************************************************

This is "butterworth.c" from the CWP/SU package with minor  
modifications by Nils Maercklin (NM, September 2008).

Reference:
Cohen, J. K. and Stockwell, J. W. (2008). CWP/SU: Seismic Un*x Release 
    No. 41: an open source software  package for seismic research and 
    processing. Center for Wave Phenomena, Colorado School of Mines. 
    http://www.cwp.mines.edu/cwpcodes/.
******************************************************************************

Copyright (c) Colorado School of Mines, 2008.
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

    *  Redistributions of source code must retain the above copyright 
       notice, this list of conditions and the following disclaimer.
    *  Redistributions in binary form must reproduce the above copyright 
       notice, this list of conditions and the following disclaimer in the 
       documentation and/or other materials provided with the distribution.
    *  Neither the name of the Colorado School of Mines nor the names of
       its contributors may be used to endorse or promote products 
       derived from this software without specific prior written permission.

Warranty Disclaimer:
THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS 
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COLORADO SCHOOL OF MINES OR CONTRIBUTORS 
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/



/*********************** self documentation **********************/
/*****************************************************************************
BUTTERWORTH - Functions to design and apply Butterworth filters:

bfdesign    design a Butterworth filter
bfhighpass  apply a high-pass Butterworth filter 
bflowpass   apply a low-pass Butterworth filter 

******************************************************************************
Function Prototypes:
void bfhighpass (int npoles, float f3db, int n, float p[], float q[]);
void bflowpass (int npoles, float f3db, int n, float p[], float q[]);
void bfdesign (float fpass, float apass, float fstop, float astop,
    int *npoles, float *f3db);

******************************************************************************
bfdesign:
Input:
fpass    frequency in pass band at which amplitude is >= apass
apass    amplitude in pass band corresponding to frequency fpass
fstop    frequency in stop band at which amplitude is <= astop
astop    amplitude in stop band corresponding to frequency fstop

Output:
npoles   number of poles
f3db     frequency at which amplitude is sqrt(0.5) (-3 db)

bfhighpass and bflowpass:
Input:
npoles   number of poles (and zeros); npoles>=0 is required
f3db     3 db frequency; nyquist = 0.5; 0.0<=f3db<=0.5 is required
n        length of p and q
p        array[n] to be filtered

Output:
q        filtered array[n] (may be equivalent to p)

******************************************************************************
Notes:
(1) Nyquist frequency equals 0.5

(2) The following conditions must be true:
    (0.0<fpass && fpass<0.5) &&
    (0.0<fstop && fstop<0.5) &&
    (fpass!=fstop) &&
    (0.0<astop && astop<apass && apass<1.0)

(3) if (fpass<fstop)

bfdesign:
Butterworth filter:  compute number of poles and -3 db frequency
for a low-pass or high-pass filter, given a frequency response
constrained at two frequencies.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/

/* #include "cwp.h" */
/* The following four lines avoid the need of "cwp.h" (NM, 09/2008). */
#include <math.h>
#ifndef PI
#define PI (3.141592653589793)
#endif

void
bfdesign (float fpass, float apass, float fstop, float astop,
    int *npoles, float *f3db)
/*****************************************************************************
Butterworth filter:  compute number of poles and -3 db frequency
for a low-pass or high-pass filter, given a frequency response
constrained at two frequencies.
******************************************************************************
Input:
fpass    frequency in pass band at which amplitude is >= apass
apass    amplitude in pass band corresponding to frequency fpass
fstop    frequency in stop band at which amplitude is <= astop
astop    amplitude in stop band corresponding to frequency fstop

Output:
npoles   number of poles
f3db     frequency at which amplitude is sqrt(0.5) (-3 db)
******************************************************************************
Notes:
(1) Nyquist frequency equals 0.5

(2) The following conditions must be true:
    (0.0<fpass && fpass<0.5) &&
    (0.0<fstop && fstop<0.5) &&
    (fpass!=fstop) &&
    (0.0<astop && astop<apass && apass<1.0)

(3) if (fpass<fstop)
        a low-pass filter is assumed
    else
        a high-pass filter is assumed
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
    float wpass,wstop,fnpoles,w3db;

    /* warp frequencies according to bilinear transform */
    wpass = 2.0*tan(PI*fpass);
    wstop = 2.0*tan(PI*fstop);

    /* if lowpass filter, then */
    if (fstop>fpass) {
        fnpoles = log((1.0/(apass*apass)-1.0)/(1.0/(astop*astop)-1.0))
            / log(pow(wpass/wstop,2.0));
        w3db = wpass/pow((1.0/(apass*apass)-1.0),0.5/fnpoles);

    /* else, if highpass filter, then */
    } else {
        fnpoles = log((1.0/(apass*apass)-1.0)/(1.0/(astop*astop)-1.0))
            / log(pow(wstop/wpass,2.0));
        w3db = wpass*pow((1.0/(apass*apass)-1.0),0.5/fnpoles);
    }

    /* determine integer number of poles */
    *npoles = 1+(int)fnpoles;

    /* determine (unwarped) -3 db frequency */
    *f3db = atan(0.5*w3db)/PI;
}

void
bfhighpass (int npoles, float f3db, int n, float p[], float q[])
/*****************************************************************************
Butterworth filter:  high-pass
******************************************************************************
Input:
npoles   number of poles (and zeros); npoles>=0 is required
f3db     3 db frequency; nyquist = 0.5; 0.0<=f3db<=0.5 is required
n        length of p and q
p        array[n] to be filtered

Output:
q        filtered array[n] (may be equivalent to p)
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
    int jpair,j;
    float r,scale,theta,a,b1,b2,pj,pjm1,pjm2,qjm1,qjm2;

    r = 2.0*tan(PI*fabs(f3db));
    if (npoles%2!=0) {
        scale = r+2.0;
        a = 2.0/scale;
        b1 = (r-2.0)/scale;
        pj = 0.0;
        qjm1 = 0.0;
        for (j=0; j<n; j++) {
            pjm1 = pj;
            pj = p[j];
            q[j] = a*(pj-pjm1)-b1*qjm1;
            qjm1 = q[j];
        }
    } else {
        for (j=0; j<n; j++)
            q[j] = p[j];
    }
    for (jpair=0; jpair<npoles/2; jpair++) {
        theta = PI*(2*jpair+1)/(2*npoles);
        scale = 4.0+4.0*r*sin(theta)+r*r;
        a = 4.0/scale;
        b1 = (2.0*r*r-8.0)/scale;
        b2 = (4.0-4.0*r*sin(theta)+r*r)/scale;
        pjm1 = 0.0;
        pj = 0.0;
        qjm2 = 0.0;
        qjm1 = 0.0;
        for (j=0; j<n; j++) {
            pjm2 = pjm1;
            pjm1 = pj;
            pj = q[j];
            q[j] = a*(pj-2.0*pjm1+pjm2)-b1*qjm1-b2*qjm2;
            qjm2 = qjm1;
            qjm1 = q[j];
        }
    }
}

void
bflowpass (int npoles, float f3db, int n, float p[], float q[])
/*****************************************************************************
Butterworth filter: low-pass
******************************************************************************
Input:
npoles   number of poles (and zeros); npoles>=0 is required
f3db     3 db frequency; nyquist = 0.5; 0.0<=f3db<=0.5 is required
n        length of p and q
p        array[n] to be filtered

Output:
q        filtered array[n] (may be equivalent to p)
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
    int jpair,j;
    float r,scale,theta,a,b1,b2,pj,pjm1,pjm2,qjm1,qjm2;

    r = 2.0*tan(PI*fabs(f3db));
    if (npoles%2!=0) {
        scale = r+2.0;
        a = r/scale;
        b1 = (r-2.0)/scale;
        pj = 0.0;
        qjm1 = 0.0;
        for (j=0; j<n; j++) {
            pjm1 = pj;
            pj = p[j];
            q[j] = a*(pj+pjm1)-b1*qjm1;
            qjm1 = q[j];
        }
    } else {
        for (j=0; j<n; j++)
            q[j] = p[j];
    }
    for (jpair=0; jpair<npoles/2; jpair++) {
        theta = PI*(2*jpair+1)/(2*npoles);
        scale = 4.0+4.0*r*sin(theta)+r*r;
        a = r*r/scale;
        b1 = (2.0*r*r-8.0)/scale;
        b2 = (4.0-4.0*r*sin(theta)+r*r)/scale;
        pjm1 = 0.0;
        pj = 0.0;
        qjm2 = 0.0;
        qjm1 = 0.0;
        for (j=0; j<n; j++) {
            pjm2 = pjm1;
            pjm1 = pj;
            pj = q[j];
            q[j] = a*(pj+2.0*pjm1+pjm2)-b1*qjm1-b2*qjm2;
            qjm2 = qjm1;
            qjm1 = q[j];
        }
    }
}
