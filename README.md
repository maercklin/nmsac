NMSAC - ANALYSIS AND PROCESSING SOFTWARE FOR SAC BINARY WAVEFORM FILES

Nils Maercklin, 2021-04-20

As part of my research work at RISSC, University of Naples Federico II 
and AMRA Scarl in Naples, Italy, I wrote this stand-alone software package 
NMSAC to read and process earthquake data in SAC format (Version 2011-07-04). 

The software is intended to be used interactively in a Unix/Linux terminal 
or within a shell script. All codes are written in C.

The codes were distributed privately throughout the years and are available 
on GitHub since April 2021.

This README contains an overview of the software, installation instructions, 
and additional information.




SOFTWARE OVERVIEW

This software collection contains programs for single- and multicomponent 
processing, for automatic onset time picking, and tools for reading and 
writing SAC header fields. The programs read and write evenly-sampled SAC 
binary waveform files (assuming NVHDR=6, IFTYPE=ITIME, LEVEN=1). Input SAC 
files are listed on the command line ("-f" option). Corresponding output 
SAC files get the same name plus an additional suffix, if the overwrite 
flag "-o" is not set. Alternatively, the programs can read from stdin and 
write to stdout. Typing any program name without further arguments prints a 
short self-documentation with parameters and defaults.

List of programs:

SACBUTTER - SAC Butterworth high-pass and/or low-pass filter
SACMAMP   - SAC Multicomponent Amplitude or RMS amplitude calculation
SACOP1    - Unary arithmetic Operations on SAC trace data
SACPOLAR  - SAC Polarization analysis of three-component data
SACPOLWH  - SAC Polarization analysis of 3-C data (Wu-Horiuchi, experimental)
SACPOFILT - SAC Polarization Filter for three-component data
SACROT3   - SAC Rotation of three-component data

SACAPKBK  - SAC Automatic single-trace P-Phase Picker (Baer-Kradolfer)
SACPKAIC  - SAC single-trace onset time estimation based on AIC
SACPKCOR  - SAC Cross-correlation based onset time pick refinement

SACLH     - List SAC Header field values 
SACSH     - Set SAC Header field values


Data processing and data analysis:

SACBUTTER is a Butterworth filter defined by the cut-off frequencies and 
the number of poles. By default, output files have the suffix ".flt", but 
input files may also be overwritten (option "-o", use with care).

SACOP1 can do various arithmetic operations on single SAC traces. These 
operations include some common mathematical functions (e.g. abs, sqrt), 
trace normalization, removal of the mean value, and approximate trace 
differentiation and integration. Output files get a suffix corresponding 
to the operation name (default), but input files may also be overwritten 
(option "-o", use with care).

SACMAMP computes the multicomponent amplitude (N-C modulus) either for 
each time sample or within a moving time window (RMS amplitude). 
Alternatively, an amplitude ratio trace can be computed.

SACPOLAR computes polarization attributes in a moving time window from 
the eigenvalues and the principal eigenvector of the three-component 
covariance matrix (e.g. Kanasewich, 1981; Jurkevics, 1988). The set of 
attributes comprises of various signal rectilinearity measures, the 
direction of polarization, and of the polarization amplitude. Several 
attributes can be computed simultaneously. Each output file has a suffix 
corresponding to the attribute name. 
SACPOLAR is somewhat similar to my program SUPOLAR (Maercklin, 1999), and 
attribute definitions can be found in the Software Manual "NMSAC-3C.pdf", 
which is also available on ResearchGate (Maercklin, 2010).

SACPOFILT is an implementation of the Montalbetti & Kanasewich (1970) 
polarization filter for three-component data. The filter utilizes a 
rectilinearity measure and the direction cosines of the principal 
polarization vector to enhance seismic phases. Before filtering the
data may be rotated into the ZRT or LQT coordinate system using 
e.g. SACROT3. Output files have the suffix ".pflt".

SACROT3 rotates three-component data into a new coordinate system, e.g. 
from ZNE into ZRT or LQT. Horizontal rotation of two-component data is 
also possible. Output files have the suffix ".rot". 

SACPOLWH is an experimental code that computes polarization attributes 
defined in Wu & Horiuchi (2008) and which may be useful to identify 
S-waves. A P-wave arrival pick must be available in the SAC header field A 
of the vertical-component trace. Output files have a suffix corresponding 
to the attribute name.


Onset time picking:

SACAPKBK is a P-phase picker for local and teleseismic events using the 
algorithm of Baer & Kradolfer (1987). A built-in Butterworth frequency 
filter may be applied before picking.

SACPKAIC assumes that the time intervals before and after the onset of a 
seismic signal are two different, locally stationary time series. The 
program analyzes the variances of the time series around an initial onset 
time estimate and searches for the minimum of the Akaike Information 
Criterion (AIC) to find the most likely separation point between the two 
time series. The time of the AIC minimum is considered to be the seismic 
phase onset (e.g. Sleeman & van Eck, 1999; Leonard & Kennett, 1999). To 
remove non-white noise from the traces, an autoregressive filter may be 
applied (prediction-error filter). Additionally, the program provides a 
frequency filter.

SACPKCOR searches for the maximum correlation coefficient between a 
reference signal (template) and the current trace in a given time window 
around an initial pick. The program may also be used as a correlation 
detector for short time segments.

SACAPKBK, SACPKAIC, SACPKCOR write the estimated onset time to the SAC 
header. The input files may be either overwritten (option "-o") or new 
files with suffix ".out" are created (default). Output waveforms are the 
unfiltered traces unless specified otherwise (option "-x").


Header utilities:

SACLH and SACSH are simple programs to read and to write header values
from and to SAC binary files. In addition to SAC waveform files, these two 
programs should also support other SAC binary file types (e.g. spectra).
SACLH writes ASCII text to stdout, whereas SACSH modifies the input SAC 
files listed after the "-f" option.




SOURCE FILES AND COMPILATION

The main programs and the libraries are written in C and listed below.


Main programs:  sacbutter.c, sacmamp.c, sacop1.c, sacrot3.c, sacpolar.c, 
                sacpolwh.c sacpofilt.c, sacpkcor.c, sacpkaic.c, sacapkbk.c, 
                saclh.c, sacsh.c
Library files:  nmsaclib.c, nmutils.c, nmpolar.c, nmetime.c, nmgeo.c
Include files:  nmsaclib.h, nmutils.h, nmsac.h, nmsacnam.h
Third-party files: butterworth.c (from CWP/SU 41).


sacbutter.c source file of SACBUTTER      (data processing/analysis)
sacmamp.c   source file of SACMAMP
sacop1.c    source file of SACOP1
sacpolar.c  source file of SACPOLAR
sacpofilt.c source file of SACPOFILT
sacrot3.c   source file of SACROT3

sacpolwh.c  source file of SACPOLWH       (experimental processing code)

sacapkbk.c  source file of SACAPKBK       (onset time picking)
sacpkaic.c  source file of SACPKAIC
sacpkcor.c  source file of SACPKCOR

saclh.c     source file of SACLH          (header utilities)
sacsh.c     source file of SACSH

nmsaclib.c  SAC-related functions used by all main programs
nmutils.c   utility functions
nmetime.c   functions for time calculations
nmgeo.c     functions for simple great-circle calculations
nmpolar.c   functions solving eigenvalue problem and for polarization analysis
butterworth.c   Functions for Butterworth low-pass/high-pass filtering 
                (taken from CWP/SU; Cohen and Stockwell, 2008).

nmsac.h     SAC binary file header definition
nmsacnam.h  global variables describing SAC header (for "nmsaclib.c" only) 
nmsaclib.h  definitions of SAC-related and other functions
nmutils.h   definitions of utility functions

Makefile    makefile for NMSAC library and main programs ("sac*.c")

To compile the main programs, edit the Makefile (if necessary) and type 
"make" in the source directory. By default, the binary executables are put 
into the "./bin" subdirectory. Assuming that the "libnmsac.a" library is 
already available you may also compile programs manually, e.g. 
"gcc -lm -lnmsac sacop1.c -o sacop1". SAC is not required to install 
and run any of these programs.




REFERENCES

Baer, M. and U. Kradolfer (1987). An automatic phase picker for local 
    and teleseismic events. Bull. Seism. Soc. Am., 77(4), 1437-1445.
Cohen, J.K. and J. W. Stockwell (2008), CWP/SU: Seismic Unix 
    Release No. 41: an open source software  package for seismic 
    research and processing, Center for Wave Phenomena, Colorado 
    School of Mines, http://www.cwp.mines.edu/cwpcodes/.
Horiuchi, S. and Y. Iio (2009). Automatic arrival time picking using
    many parameters for the onset discrimination.
    Second Earthquake Early Warning Workshop, Tokyo, April 21-22.
Jurkevics, A. (1988). Polarization analysis of three-component array data.
    Bull. Seism. Soc. Am., 78(5), 1725-1743.
Kanasewich, E. R. (1981). Time Sequence Analysis in Geophysics.
    The University of Alberta Press, ISBN 0888640749.
Leonard, M. and B. L. N. Kennett (1999). Multi-component autoregressive
    techniques for the analysis of seismograms. Phys. Earth Planet.
    Int., 113, 247-263, doi:10.1016/S0031-9201(99)00054-0.
Maercklin, N. (1999). Polarisationsanalyse refraktionsseismischer Daten vom
    Vulkan Merapi, Indonesien. Masters thesis (Diplomarbeit), Geophysics,
    Christian-Albrechts-University, Kiel, Germany.
Maercklin, N. (2010). Three-component processing and analysis tools for seismic 
    data in SAC format. Software manual available on ResearchGate, 
    doi:10.13140/2.1.2222.0801.
Montalbetti, J. R. and E. R. Kanasewich (1970). Enhancement of teleseismic
    body phases with a polarization filter. Geophys. J. R. Astr. Soc., 
    21(2), 119-129, doi:10.1111/j.1365-246X.1970.tb01771.x.
Sleeman, R. and T. van Eck (1999). Robust automatic P-phase picking:
    an on-line implementation in the analysis of broadband seismogram
    recordings. Phys. Earth Planet. Int., 113, 265-275, 
    doi:10.1016/S0031-9201(99)00007-2.
Wu, C. and S. Horiuchi (2008). Automatic determination of source parameters
    of the 2007 Noto Hanto earthquake. Earth, Planets, and Space, 
    60, 1053-1057.

The SAC (Seismic Analysis Code) user's manual and further information 
on SAC is available at http://www.iris.edu/manuals/sac/manual.html.




LICENSE

Copyright (c) 2021 Nils Maercklin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.




CONTACT

Nils Maercklin
formerly at: RISSC, University of Naples Federico II, Naples, Italy
https://www.researchgate.net/profile/Nils-Maercklin
https://www.linkedin.com/in/maercklin/
