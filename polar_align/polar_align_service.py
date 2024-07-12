from device import mount, ccd
import PyIndi


from astropy.coordinates import FK5, AltAz
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.wcs import WCS
from astropy.io import fits
from datetime import datetime
from functools import lru_cache
import numpy as np
import os.path as op
import scipy.optimize
import subprocess
import traceback

class PolarAlignService:

    indiclient: PyIndi.BaseClient
    telescope: mount.Mount
    camera:  ccd.Ccd
    lat, lon = 50.1, 14.4
    site = EarthLocation.from_geodetic(lon=lon, lat=lat)
    image_width_pix = 1280
    image_height_pix = 1024
    first_capture = ""
    next_capture = ""
    NCP_to_date_altaz_deg = []
    R0_altaz_deg = []
    polar_align_started = False

    def __init__(self, mount, ccd, indiclient) -> None:
        super().__init__()
        self.telescope = mount
        self.camera = ccd
        self.indiclient = indiclient

    def StartPolarAlign(self, request, context):
        print("called")
        try:
            print("start polar align")
            self.polar_align_started = True
            self.camera.setExposure(2)
            self.first_capture = self.camera.takeExposure()
            self.telescope.goto("3h", 0)
        
            # step (1)
            A_J2000 = self.get_J2000(self.first_capture)
            A_pix = self.J2000_to_pix(self.first_capture, A_J2000)
            A_time = self.get_time(self.first_capture)
            A_altaz_frame = self.get_altaz_frame(self.first_capture)

            self.camera.setExposure(2)
            self.next_capture = self.camera.takeExposure()

            # step (2)
            B_J2000 = self.get_J2000(self.next_capture)
            #B_pix = self.J2000_to_pix(self.first_capture, B_J2000)
            #B_altaz_frame = self.get_altaz_frame(self.next_capture)
            #B_altaz = self.get_altaz(self.next_capture)
            #B_altaz_deg = self.get_altaz_deg(self.next_capture)
            B_time = self.get_time(self.next_capture)

            # step (3)
            dt_BA_sec = (B_time - A_time).total_seconds()

            
            res = scipy.optimize.minimize(self.displacement_error_squared_time_comp, x0=[0,0], args=(A_pix, A_time, dt_BA_sec), method="Nelder-Mead")
            print("res ", res)
            print("res x ", res.x)
            
            R0_pix = res.x
            R0_J2000 = self.pix_to_J2000(self.first_capture, R0_pix)
            R0_altaz = R0_J2000.transform_to(A_altaz_frame)
            self.R0_altaz_deg = np.asarray([R0_altaz.alt.deg, R0_altaz.az.deg])
            print("R_0:", R0_pix, R0_J2000, sep="\n", end="\n\n")

            NCP_J2000 = SkyCoord(ra=0, dec=90, frame='fk5', unit="deg", equinox="J2000")
            NCP_to_date_J2000 = SkyCoord(ra=0, dec=90, frame='fk5', unit="deg", equinox=self.get_time(self.first_capture)).transform_to(
                FK5(equinox="J2000"))
            NCP_to_date_pix = self.J2000_to_pix(self.first_capture, NCP_to_date_J2000)
            NCP_to_date_altaz = NCP_to_date_J2000.transform_to(A_altaz_frame)
            self.NCP_to_date_altaz_deg = np.asarray([NCP_to_date_altaz.alt.deg, NCP_to_date_altaz.az.deg])
            #return super().StartPolarAlign(request, context)
        except Exception as err:
            print("failed to polar align: ", err)
            print(traceback.print_exc())
            self.telescope.Park()
            exit(0)

    def StopPolarAlign(self, request, context):
        return super().StopPolarAlign(request, context)
    
    def ComputePolarError(self):
        if not(self.polar_align_started):
            raise Exception("Polar align has not started")
        print("compute polar error")

        # steps (5), (6)
        self.camera.setExposure(2)
        I=self.camera.takeExposure()
        I_time = self.get_time(I)
        I_J2000 = self.get_J2000(I)
        I_altaz = self.get_altaz(I)
        I_altaz_deg = self.get_altaz_deg(I)
        #output_Ii_altaz_deg.append(I_altaz_deg)

        offset_altaz_deg = I_altaz_deg - self.get_altaz_deg(self.next_capture)
        Ri_altaz_deg = self.R0_altaz_deg + offset_altaz_deg
        # output_Ri_altaz_deg.append(Ri_altaz_deg)

        p1 = SkyCoord(Ri_altaz_deg[1], Ri_altaz_deg[0], unit="deg")
        p2 = SkyCoord(self.NCP_to_date_altaz_deg[1], self.NCP_to_date_altaz_deg[0], unit="deg")
        error_deg = p1.separation(p2).deg
        error_az, error_alt = p2.spherical_offsets_to(p1)

        error_arcmin = error_deg * 60
        # output_x_errors_arcmin.append(error_az.deg * 60)
        #output_y_errors_arcmin.append(error_alt.deg * 60)
        #output_tot_errors_arcmin.append(error_arcmin)
        print("azimuth error ", error_az.deg , " alt error ", error_alt.deg )


    def get_wcs_header(self, name: str):
        if not(op.exists(f"{name}.wcs")):
            print(f"platesolving image {name}.fits")
            subprocess.run(["astap_cli", "-f",  f"{name}.fits", "-r", "15", "-fov", "0.5", "-wcs", "-speed", "slow"], check=True, text=True)
        path = f"{name}.wcs"
        with fits.open(path, "readonly", encoding="ascii", errors="ignore") as fp:
            #return fits.Header.fromstring(fp.read(), sep='\n')
            return fp[0].header

    @lru_cache(maxsize=None)
    def get_wcs(self, name: str):
        return WCS(self.get_wcs_header(name))

    def get_time(self, name: str) -> datetime:
        h = self.get_wcs_header(name)
        timestamp_str = h["DATE-OBS"]
        return datetime.fromisoformat(timestamp_str)

    def get_J2000(self, name: str) -> SkyCoord:
        h = self.get_wcs_header(name)
        ra, dec = h["RA"], h["DEC"]
        coord = SkyCoord(ra=ra, dec=dec, frame='fk5', unit="deg", equinox="J2000")
        return coord

    def get_altaz(self, name: str) -> SkyCoord:
        coord = self.get_J2000(name)
        time = self.get_time(name)
        altaz_frame = AltAz(obstime=time, location=self.site)
        return coord.transform_to(altaz_frame)

    def get_altaz_deg(self, name: str) -> np.ndarray:
        coord = self.get_altaz(name)
        return np.asarray([coord.alt.deg, coord.az.deg])

    def get_altaz_frame(self, name: str) -> AltAz:
        coord = self.get_J2000(name)
        time = self.get_time(name)
        altaz_frame = AltAz(obstime=time, location=self.site)
        return altaz_frame

    def J2000_to_pix(self, reference_name: str, coord: SkyCoord) -> np.ndarray:
        wcs = self.get_wcs(reference_name)
        return np.asarray(wcs.world_to_pixel(coord))

    def pix_to_J2000(self, reference_name: str, pix: np.ndarray) -> SkyCoord:
        x, y = pix
        wcs = self.get_wcs(reference_name)
        return wcs.pixel_to_world(x, y)

    def J2000_advance_time(self, coord: SkyCoord, dt_sec: float) -> SkyCoord:
        return SkyCoord(ra=coord.ra.deg - 360/24/60/60 * dt_sec,
                        dec=coord.dec.deg, frame='fk5', unit="deg", equinox="J2000")

    def displacement_error_squared_time_comp(self, x, pix, A_time, dt_BA_sec):
        coord_J2000_A = self.pix_to_J2000(self.first_capture, pix)
        
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
            
        pix2 = self.J2000_to_pix(self.next_capture, coord_J2000_B)
        return np.sum((pix - pix2) ** 2)

