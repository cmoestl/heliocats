Numpy record arrays of solar wind in situ data in python pickle

authors: C. Möstl, Austrian Space Weather Office, Graz, Austria (twitter @chrisoutofspace, chris.moestl@outlook.com), A. J. Weiss, M. A. Reiss, R. L. Bailey,  A. Isavnin, R. M. Winslow, C. L. Simon Wedlund, D. Stansby, A. Isavnin
made with https://github.com/cmoestl/heliocats

files available:

active spacecraft (status April 2023):

Solar Orbiter

only magnetic field:
solo_2020_april_2022_sep_mag_sceq.p
solo_2020_april_2022_sep_mag_rtn.p


Parker Solar Probe until 2021 Dec 30:

psp_2018_2022_rtn.p
psp_2018_2022_sceq.p

only magnetic field from 2021 Dec 31 onwards:
(this file has no header)
psp_2022_add_mag_rtn.p
psp_2022_add_mag_sceq.p



STEREO-A science data
stereoa_2007_2020_sceq.p
stereoa_2007_2020_rtn.p
Beacon data
stereoa_2020_now_rtn_beacon.p
stereoa_2020_now_sceq_beacon.p

Wind:
wind_1995_2021_heeq.p
wind_1995_2021_gse.p

more recent data are included in:
wind_2018_now_heeq.p
wind_2018_now_gse.p

MAVEN
maven_2014_2018_removed_smoothed.p
maven_2014_2018_removed.p
maven_2014_2018.p


finished missions:

ulysses_1990_2009_rtn.p

stereob_2007_2014_sceq.p
stereob_2007_2014_rtn.p

messenger_2007_2015_sceq_removed.p
messenger_2007_2015_sceq.p

vex_2007_2014_sceq_removed.p
vex_2007_2014_sceq.p



Load into python with e.g.:

>> import pickle
>> file='stereob_2007_2014_sceq.p'
>> [data, header]=pickle.load(open(file,'rb'))


see header by
>> print(header)


access data by e.g.
>> data.time
>> data.bt
>> data.vt
as "data" is a numpy recarray


variables are:
time	    datetime object	
bt          total magnetic field
bxyz        magnetic field components
vt          total proton speed
np          proton density
tp          proton temperature
xyz         sapcecraft position in HEEQ cartesian coordinates
r/lat/lon   spacecraft position in HEEQ spherical coordinates

For units see header string.

Spacecraft positions are given in Heliocentric Earth Equatorial Coordinates (HEEQ) coordinates.
Here, Z is the solar rotation axis, X the intersection of the solar equator and solar central 
meridian as seen from Earth, and Y completes the right handed triad.
 
See Hapgood 1992 http://www.igpp.ucla.edu/public/vassilis/ESS261/Lecture03/Hapgood_sdarticle.pdf 
for how to make conversions between coordinate systems.

Coordinate system for all magnetic field components is SCEQ except for Wind (HEEQ), Ulysses (RTN) and MAVEN (MSO).
Definition of SpaceCraft Equatorial Coordinates (SCEQ):
      Z is the solar rotation axis.
      Y is the cross product of Z and R, with R being the vector that points from the Sun to the spacecraft.
      X completes the right handed triad.
      This system is thus like HEEQ but centered on the respective in situ spacecraft, so 
      the SCEQ X and Y base vectors are rotated by the HEEQ longitude of the in situ spacecraft from HEEQ X and Y.
      The Y vector is similar to the T vector in an RTN system for each spacecraft, but the X and Z vectors 
      are rotated around Y compared to an RTN system.
      We choose SCEQ because it has the advantage that a comparison between multipoint CME events 
      and for comparison to simulations there is always a similar reference plane (the solar equatorial plane).


Example of loaded header string and data array for file='stereob_2007_2014_sceq_ndarray.p':

header:

'STEREO-B magnetic field (IMPACT instrument) and plasma data (PLASTIC, science), obtained from https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level2/behind/magplasma   Timerange: 2007-Jan-01 00:00 to 2014-Sep-27 16:31, with a an average time resolution of 61 seconds. The data are available in a numpy recarray, fields can be accessed by stb.time, stb.bx, stb.vt etc. Missing data has been set to "np.nan". Total number of data points: 4001229. Units are btxyz [nT, SCEQ], vt [km/s], np[cm^-3], tp [K], heliospheric position x/y/z/r/lon/lat [AU, degree, HEEQ]. Made with https://github.com/cmoestl/heliocats and https://github.com/heliopython/heliopy. By C. Moestl (twitter @chrisoutofspace), A. J. Weiss, R. L. Bailey and D. Stansby. File creation date: 2020-Apr-06 16:54 UTC'


data:

rec.array([(datetime.datetime(2007, 1, 1, 0, 0), -0.17943429, -2.18674636,  1.7932473 , 2.83691227, nan, nan, nan,  0.97666821,  0.00225768, -0.05210046, 0.97805948, -3.05354637,  1.32445475e-01),
           (datetime.datetime(2007, 1, 1, 0, 1), -0.13964184, -2.18558764,  1.65750451, 2.7494069 , nan, nan, nan,  0.97666808,  0.00225762, -0.05210185, 0.97805943, -3.05362816,  1.32442261e-01),
           (datetime.datetime(2007, 1, 1, 0, 2), -0.16222884, -2.02721858,  1.57420026, 2.57451478, nan, nan, nan,  0.97666796,  0.00225757, -0.05210324, 0.97805938, -3.05370994,  1.32439047e-01),
           ...,
           (datetime.datetime(2014, 9, 27, 16, 29), -3.03044344,  0.2601788 , -1.07245173, 3.22964609, nan, nan, nan, -0.97304608, -0.33086124, -0.10771433, 1.03338773, -5.98304654, -1.61220617e+02),
           (datetime.datetime(2014, 9, 27, 16, 30), -3.09114393,  0.49015975, -1.21775023, 3.36384756, nan, nan, nan, -0.97304679, -0.33086091, -0.10771352, 1.03338821, -5.98299865, -1.61220647e+02),
           (datetime.datetime(2014, 9, 27, 16, 31),         nan,         nan,         nan,        nan, nan, nan, nan, -0.9730475 , -0.33086058, -0.10771271, 1.03338869, -5.98295076, -1.61220677e+02)],
          dtype=[('time', 'O'), ('bx', '<f8'), ('by', '<f8'), ('bz', '<f8'), ('bt', '<f8'), ('vt', '<f8'), ('np', '<f8'), ('tp', '<f8'), ('x', '<f8'), ('y', '<f8'), ('z', '<f8'), ('r', '<f8'), ('lat', '<f8'), ('lon', '<f8')])


