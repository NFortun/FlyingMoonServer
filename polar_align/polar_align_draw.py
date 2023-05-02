from astropy.coordinates import FK5, AltAz
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.wcs import WCS
from astropy.io import fits
from datetime import datetime
from functools import lru_cache
from matplotlib import ticker, cm
import configparser
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os.path as op
import pandas as pd
import scipy.optimize
import subprocess

lat, lon = 50.1, 14.4
site = EarthLocation.from_geodetic(lon=lon, lat=lat)
arcsec_per_pix = 2.01  # TODO it's best not to rely on this, it's approximate
image_width_pix = 1936
image_height_pix = 1088

def get_image(name: str):
    path = f"./01-captured-images/{name}.png"
    return plt.imread(path)

def get_wcs_header(name: str):
    path = f"./03-astap-analysis/{name}.wcs"
    with open(path, "r", encoding="ascii", errors="ignore") as fp:
        return fits.Header.fromstring(fp.read(), sep='\n')

@lru_cache(maxsize=None)
def get_wcs(name: str):
    return WCS(get_wcs_header(name))

def get_time(name: str) -> datetime:
    path = f"./01-captured-images/{name}.CameraSettings.txt"
    parser = configparser.ConfigParser()
    parser.read(path)
    timestamp_str = parser["ZWO ASI290MM Mini"]["TimeStamp"][:26]
    return datetime.fromisoformat(timestamp_str)

def get_J2000(name: str) -> SkyCoord:
    h = get_wcs_header(name)
    ra, dec = h["CRVAL1"], h["CRVAL2"]
    coord = SkyCoord(ra=ra, dec=dec, frame='fk5', unit="deg", equinox="J2000")
    return coord

def get_altaz(name: str) -> SkyCoord:
    coord = get_J2000(name)
    time = get_time(name)
    altaz_frame = AltAz(obstime=time, location=site)
    return coord.transform_to(altaz_frame)

def get_altaz_deg(name: str) -> np.ndarray:
    coord = get_altaz(name)
    return np.asarray([coord.alt.deg, coord.az.deg])

def get_altaz_frame(name: str) -> AltAz:
    coord = get_J2000(name)
    time = get_time(name)
    altaz_frame = AltAz(obstime=time, location=site)
    return altaz_frame

def J2000_to_pix(reference_name: str, coord: SkyCoord) -> np.ndarray:
    wcs = get_wcs(reference_name)
    return np.asarray(wcs.world_to_pixel(coord))

def pix_to_J2000(reference_name: str, pix: np.ndarray) -> SkyCoord:
    x, y = pix
    wcs = get_wcs(reference_name)
    return wcs.pixel_to_world(x, y)

def J2000_advance_time(coord: SkyCoord, dt_sec: float) -> SkyCoord:
    return SkyCoord(ra=coord.ra.deg - 360/24/60/60 * dt_sec,
                    dec=coord.dec.deg, frame='fk5', unit="deg", equinox="J2000")

def get_image_corners_J2000(name: str):
    return [
        pix_to_J2000(name, (0, 0)),
        pix_to_J2000(name, (image_width_pix, 0)),
        pix_to_J2000(name, (image_width_pix, image_height_pix)),
        pix_to_J2000(name, (0, image_height_pix)),
    ]

all_image_names = [f"_{i:05}" for i in range(3, 26)]
A = all_image_names[0]
B = all_image_names[1]



# All values "*_pix" are given in A local coordinates (ie. pixels)
output_Ii_altaz_deg = []
output_Ri_altaz_deg = []
output_x_errors_arcmin = []
output_y_errors_arcmin = []
output_tot_errors_arcmin = []
output_xy_is_altaz = True

compensate_AB_time = True  # <--- this should be True for best results

# step (1)
A_J2000 = get_J2000(A)
A_pix = J2000_to_pix(A, A_J2000)
A_time = get_time(A)
A_altaz_frame = get_altaz_frame(A)

# step (2)
B_J2000 = get_J2000(B)
B_pix = J2000_to_pix(A, B_J2000)
B_altaz_frame = get_altaz_frame(B)
B_altaz = get_altaz(B)
B_altaz_deg = get_altaz_deg(B)
B_time = get_time(B)

# step (3)
dt_BA_sec = (B_time - A_time).total_seconds()

def displacement_error_squared_time_comp(pix):
    coord_J2000_A = pix_to_J2000(A, pix)
    
    if compensate_AB_time:
        # What we really want here is to project A(pix) -> Alt/Az -> B(pix2),
        # but platesolving gives us the projection in terms of Ra/Dec J2000.
        # Simply round-tripping through J2000 is not quite correct, since
        # A and B are taken at different times. (Imagine if the mount didn't
        # move at all between A and B - they would have same Alt/Az coordinates,
        # but sightly different J2000.) What we can do is to "fix" the J2000
        # coordinates, baking in the rotation due to time before looking them up in B.
        # Note that we rotate (ie. increment RA) around the CP to date, not the J2000 one.
        coord_to_date_A = coord_J2000_A.transform_to(FK5(equinox=A_time))
        coord_to_date_B = SkyCoord(ra=coord_to_date_A.ra.deg + 360 / 24 / 60 / 60 * dt_BA_sec,
                                   dec=coord_to_date_A.dec.deg,  # ^ approx. 360°/day
                                   frame='fk5', unit="deg", equinox=A_time)
        coord_J2000_B = coord_to_date_B.transform_to(FK5(equinox="J2000"))
    else:
        coord_J2000_B = coord_J2000_A
        
    pix2 = J2000_to_pix(B, coord_J2000_B)
    return np.sum((pix - pix2) ** 2)


res = scipy.optimize.minimize(displacement_error_squared_time_comp, A_pix, method="nelder-mead")
R0_pix = res.x
R0_J2000 = pix_to_J2000(A, R0_pix)
R0_altaz = R0_J2000.transform_to(A_altaz_frame)
R0_altaz_deg = np.asarray([R0_altaz.alt.deg, R0_altaz.az.deg])
print("R_0:", R0_pix, R0_J2000, sep="\n", end="\n\n")

# step (4)
NCP_J2000 = SkyCoord(ra=0, dec=90, frame='fk5', unit="deg", equinox="J2000")
NCP_to_date_J2000 = SkyCoord(ra=0, dec=90, frame='fk5', unit="deg", equinox=get_time(A)).transform_to(
    FK5(equinox="J2000"))
NCP_to_date_pix = J2000_to_pix(A, NCP_to_date_J2000)
NCP_to_date_altaz = NCP_to_date_J2000.transform_to(A_altaz_frame)
NCP_to_date_altaz_deg = np.asarray([NCP_to_date_altaz.alt.deg, NCP_to_date_altaz.az.deg])

# steps (5), (6)
for I in all_image_names[1:]:
    I_time = get_time(I)
    I_J2000 = get_J2000(I)
    I_altaz = get_altaz(I)
    I_altaz_deg = get_altaz_deg(I)
    output_Ii_altaz_deg.append(I_altaz_deg)

    offset_altaz_deg = I_altaz_deg - B_altaz_deg
    Ri_altaz_deg = R0_altaz_deg + offset_altaz_deg
    output_Ri_altaz_deg.append(Ri_altaz_deg)

    p1 = SkyCoord(Ri_altaz_deg[1], Ri_altaz_deg[0], unit="deg")
    p2 = SkyCoord(NCP_to_date_altaz_deg[1], NCP_to_date_altaz_deg[0], unit="deg")
    error_deg = p1.separation(p2).deg
    error_az, error_alt = p2.spherical_offsets_to(p1)

    error_arcmin = error_deg * 60
    output_x_errors_arcmin.append(error_az.deg * 60)
    output_y_errors_arcmin.append(error_alt.deg * 60)
    output_tot_errors_arcmin.append(error_arcmin)

print("Final error [arcmin]:", output_tot_errors_arcmin[-1])

pd.DataFrame.from_dict(dict(i=all_image_names[1:],
                            Ii_altaz_deg=output_Ii_altaz_deg,
                            Ri_altaz_deg=output_Ri_altaz_deg,
                            tot_error_arcmin=output_tot_errors_arcmin)).set_index("i")

