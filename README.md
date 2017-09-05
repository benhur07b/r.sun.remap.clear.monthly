### REMap-SOLAR monthly average clear-sky GHI script

MODULE: 
r.sun.remap.clear.monthly.py

AUTHOR(S):
Ben Hur S. Pintor bhs.pintor@gmail.com

PURPOSE:
This script computes for the monthly clear sky GHI values using an elevation raster and monthly average linke turbidity coefficients.
The elevation raster is used to compute for the slope, aspect, and horizon maps.
The monthly linke turbidity coefficients is a .csv file
The monthly average is taken as the GHI value during the average day of the month (Duffie and Beckman, 2013).

COPYRIGHT:    (C) 2016 by Ben Hur S. Pintor

This program is free software under the GNU General Public License (>=v2). Read the file COPYING that comes with GRASS for details and the LICENSE that comes with this file.