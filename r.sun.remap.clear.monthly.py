#!/usr/bin/env python
#
##############################################################################
#
# MODULE:   r.sun.remap.clear.monthly.py
# AUTHOR(S):    Ben Hur S. Pintor <bhs dot pintor at gmail dot com>
#
# PURPOSE:  This script computes for the monthly clear sky GHI values using
#       an elevation raster and monthly average linke turbidity
#       coefficients.
#
#       The elevation raster is used to compute for the slope, aspect, and
#       horizon maps.
#
#       The monthly linke turbidity coefficients is a .csv file
#
#       The monthly average is taken as the GHI value during the average
#       day of the montth (Duffie and Beckman, 2013).
#
#
# COPYRIGHT:    (C) 2016 by Ben Hur S. Pintor
#
#       This program is free software under the GNU General Public
#       License (>=v2). Read the file COPYING that comes with GRASS
#       for details.
#
##############################################################################
#
#%Module
#% description: Computes for the monthly clear-sky GHI using an elevation raster and monthly average linke turbidity coefficient .csv (For REMap SOLAR Component)
#% keyword: raster
#% keyword: solar radiation
#%end
#%flag
#% key: h
#% description: Remove horizon rasters when the process ends
#%end
#%option
#% key: elevation
#% type: string
#% description: An elevation raster file
#% gisprompt: new,bin,file
#% required: yes
#%end
#%option
#% key: linke
#% type: string
#% description: A .csv file of monthly linke turbidity coefficients
#% gisprompt: new,csv,file
#% required: yes
#%end
#%option
#% key: prefix
#% type: string
#% description: The prefix to append to the maps (i.e. province name, etc)
#% required: yes
#%end
#%option
#% key: step
#% type: double
#% description: Angle step size for horizon computations [degrees]
#% required : no
#% answer: 3.75
#%end
#%option
#% key: bufferzone
#% type: double
#% description: The bufferzone around present region for the horizon computations
#% required : no
#% answer: 10000.0
#%end
#%option
#% key: maxdistance
#% type: double
#% description: The maximum distance for the horizon computations
#% required : no
#% answer: 10000.0
#%end
#%option
#% key: timestep
#% type: double
#% description: The time step for computing solar radiation using r.sun
#% required : no
#% answer: 0.08333
#%end
#%option
#% key: npartitions
#% type: integer
#% description: The number of chunks to divide the DEM while running r.sun
#% required : no
#% answer: 4
#%end
#%option
#% key: tension
#% type: double
#% description: The tension parameter for linke interpolation
#% required : no
#% answer: 100.0
#%end
#%option
#% key: smooth
#% type: double
#% description: The smoothing for linke interpolation
#% required : no
#% answer: 0.0
#%end
#%option
#% key: resolution
#% type: integer
#% description: The horizontal resolution for the horizon computations
#% required : no
#% answer: 100
#%end



import os
import sys
import grass.script as gscript
from grass.pygrass.gis import Mapset
from grass.pygrass.modules import Module


'''Mean days list (from Duffie and Beckman, 1991)'''
MEANDAYS = [{'month': 'JAN', 'julian day': 17, 'declination': -0.364773813667},
            {'month': 'FEB', 'julian day': 47, 'declination': -0.226892802759},
            {'month': 'MAR', 'julian day': 75, 'declination': -0.0418879020479},
            {'month': 'APR', 'julian day': 105, 'declination': 0.164060949687},
            {'month': 'MAY', 'julian day': 135, 'declination': 0.328121899375},
            {'month': 'JUN', 'julian day': 162, 'declination': 0.403171057211},
            {'month': 'JUL', 'julian day': 198, 'declination': 0.370009801423},
            {'month': 'AUG', 'julian day': 228, 'declination': 0.235619449019},
            {'month': 'SEP', 'julian day': 258, 'declination': 0.0383972435439},
            {'month': 'OCT', 'julian day': 288, 'declination': -0.167551608191},
            {'month': 'NOV', 'julian day': 318, 'declination': -0.329867228627},
            {'month': 'DEC', 'julian day': 344, 'declination': -0.401425727959}]


'''Create objects for the GRASS modules to be used'''
r_in_gdal = Module("r.in.gdal")
r_slope_aspect = Module("r.slope.aspect")
r_horizon = Module("r.horizon")
r_mapcalc = Module("r.mapcalc")
r_sun = Module("r.sun")
r_resamp_interp = Module("r.resamp.interp")

v_in_ascii = Module("v.in.ascii")
v_db_addcolumn = Module("v.db.addcolumn")
v_what_rast = Module("v.what.rast")
v_db_update = Module("v.db.update")
v_surf_rst = Module("v.surf.rst")
v_in_ogr = Module("v.in.ogr")

g_region = Module("g.region")
g_remove = Module("g.remove")

m_proj = Module("m.proj")


def in_mapset(m, dtype, pattern=None):
    """Checks if a map with the same name already exists in the mapset"""

    mset = Mapset()
    maps = mset.glist(dtype, pattern=pattern)

    if m in maps:
        return True
    else:
        return False


def cleanup_horizon(prefix):
    """Removes the horizon rasters when processing is over. Done to save space."""

    g_remove(flags="rf", type="raster", pattern="{}_horizon".format(prefix))


def import_elev(elevation, prefix):

    if in_mapset("{}_elev".format(prefix), "raster"):
        print "{}_elev raster already exists".format(prefix)
        pass

    else:
        r_in_gdal(input=elevation, output="{}_elev".format(prefix))


def slope_and_aspect(prefix):

    if in_mapset("{}_slope".format(prefix), "raster") and in_mapset("{}_aspect".format(prefix), "raster"):
        print "{}_slope and {}_aspect rasters already exist".format(prefix, prefix)
        pass

    else:
        r_slope_aspect(elevation="{}_elev".format(prefix), slope="{}_slope".format(prefix), aspect="{}_aspect".format(prefix), overwrite=True)


def prepare_linke(linke, prefix):

    linkeproj = "{}_proj.csv".format(linke[:-4])
    m_proj(input=linke, output=linkeproj, separator="comma", flags="i", overwrite=True)

    if in_mapset("{}_linke".format(prefix), "vector"):
        print "{}_linke vector already exists".format(prefix)
        pass

    else:
        mtxt = ""
        for m in MEANDAYS:
            month = m['month']
            mtxt += ",{} DOUBLE PRECISION".format(month)

        cols = "X DOUBLE PRECISION,Y DOUBLE PRECISION{}".format(mtxt)

        v_in_ascii(input=linkeproj, output="{}_linke".format(prefix), separator=", space",
            columns=cols)

        # Add columns for the normalized Linke values
        mtxt = ""
        for m in MEANDAYS:
            month = m['month']
            mtxt += ",{}_NORM DOUBLE PRECISION".format(month)

        cols = "ELEVATION DOUBLE PRECISION{}".format(mtxt)

        v_db_addcolumn(map="{}_linke".format(prefix), columns=cols)

        # Sample the elevation raster at the Linke points
        v_what_rast(map="{}_linke".format(prefix), raster="{}_elev".format(prefix), column="ELEVATION")

        # Compute for the normalized Linke turbidity values
        for m in MEANDAYS:
            month = m['month']
            col = "{}_NORM".format(month)
            update = "{} + (0.00035 * {})".format(month, "ELEVATION")
            v_db_update(map="{}_linke".format(prefix), layer=1, column=col, query_column=update)


def normalized_linke_rasters(tension, smooth, prefix):
    """Computes for the monthly normalized Linke turbidity rasters (12)."""

    for m in MEANDAYS:

        month = m['month']
        col = "{}_NORM".format(month)
        name = "{}_linke_{}".format(prefix,month)

        if in_mapset(name, "raster"):
            print "{} already exists".format(name)
            pass

        else:
            v_surf_rst(input="{}_linke".format(prefix), zcolumn=col, elevation=name, mask="{}_elev".format(prefix),
                tension=tension, smooth=smooth, overwrite=True)

            print "Done with Normalized Linke for {}".format(m)


def prepare_horizon(elevation, step, bufferzone, maxdistance, resolution, prefix):

    g_region(raster="{}_elev".format(prefix), res=resolution)

    mset = Mapset()
    maps = mset.glist("raster")
    hors = [h for h in maps if "{}_horizon".format(prefix) in h]

    if len(hors) == int(360/float(step)):
        pass

    else:
        r_horizon(elevation="{}_elev".format(prefix), step=step, bufferzone=bufferzone, maxdistance=maxdistance, output="{}_horizon".format(prefix), overwrite=True)


def prepare_rsun_inputs(elevation, linke, step, bufferzone, maxdistance, tension, smooth, resolution, prefix):

    import_elev(elevation, prefix)

    g_region(raster="{}_elev".format(prefix), res=10)
    slope_and_aspect(prefix)

    prepare_linke(linke, prefix)
    normalized_linke_rasters(tension, smooth, prefix)

    prepare_horizon(elevation, step, bufferzone, maxdistance, resolution, prefix)


def monthly_clear_sky_ghi(elevation, linke, step, bufferzone, maxdistance, tension, smooth, resolution, timestep, npartitions, prefix):

    prepare_rsun_inputs(elevation, linke, step, bufferzone, maxdistance, tension, smooth, resolution, prefix)

    g_region(raster="{}_elev".format(prefix), res=10)

    for m in MEANDAYS:
        month = m["month"]
        date = m["julian day"]
        decl = m["declination"]

        if in_mapset("{}_GHI_clear_{}".format(prefix, month), "raster"):
            print "{}_GHI_clear_{} already in mapset".format(prefix, month)
            pass

        else:
            r_sun(elevation="{}_elev".format(prefix),
                  aspect="{}_aspect".format(prefix),
                  slope="{}_slope".format(prefix),
                  linke="{}_linke_{}".format(prefix, month),
                  horizon_basename="{}_horizon".format(prefix),
                  horizon_step=step,
                  glob_rad="{}_GHI_clear_{}".format(prefix, month),
                  day=date,
                  step=timestep,
                  declination=decl,
                  npartitions=npartitions)



def main():

    '''Set options'''
    elevation = options['elevation']
    linke = options['linke']
    step = options['step']
    bufferzone = options['bufferzone']
    maxdistance = options['maxdistance']
    timestep = options['timestep']
    npartitions = options['npartitions']
    tension = options['tension']
    smooth = options['smooth']
    resolution = options['resolution']
    prefix = options['prefix']

    '''Set flags'''
    flag_h = flags['h']

    monthly_clear_sky_ghi(elevation, linke, step, bufferzone, maxdistance, tension, smooth, resolution, timestep, npartitions, prefix)

    if flag_h:
        cleanup_horizon(prefix)

    return 0

if __name__ == "__main__":
    options, flags = gscript.parser()
    # atexit.register(cleanup_horizon)
    sys.exit(main())
