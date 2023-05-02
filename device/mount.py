import time
import PyIndi

from astropy.coordinates import SkyCoord
from astropy import units as u
from device.base_device import BaseDevice

class Mount(BaseDevice):
    def __init__(self, name, mount_device, indiclient) -> None:
        super().__init__(name, mount_device, indiclient)

    def goto(self, ra=0, dec=0):
        on_coord_set=self.getSwitch("ON_COORD_SET")
       # on_coord_set[0].s=PyIndi.ISS_OFF  # TRACK
       # on_coord_set[1].s=PyIndi.ISS_ON # SLEW
       # on_coord_set[2].s=PyIndi.ISS_OFF # SYNC
        on_coord_set.setRule("TRACK")
        on_coord_set.setState(PyIndi.ISS_ON)
        self.indiclient.sendNewSwitch(on_coord_set)

        equatorial_eod_coord = self.getNumber("EQUATORIAL_EOD_COORD")
        coordinates = SkyCoord(ra=ra, dec=dec, frame="icrs", unit=(u.hour, u.deg))
        if(ra):
            equatorial_eod_coord[0].setValue(coordinates.ra.hour)
        if(dec):
            equatorial_eod_coord[1].setValue(coordinates.dec.value)
        self.indiclient.sendNewNumber(equatorial_eod_coord)
        while(equatorial_eod_coord.getState()==PyIndi.IPS_BUSY):
            time.sleep(2)
            print("Telescope is moving: ", equatorial_eod_coord[0].value, " ", equatorial_eod_coord[1].value)
        print("Settling...")
        time.sleep(5)

    def Park(self):
        print("parking mount")
        park=self.getSwitch("TELESCOPE_PARK")
        while(not(park)):
            park=self.getSwitch("TELESCOPE_PARK")
        park[0].setState(PyIndi.ISS_ON)
        park[1].setState(PyIndi.ISS_OFF)
        self.indiclient.sendNewSwitch(park)

    def unpark(self):
        # probleme
        print("unparking mount")
        unpark=self.getSwitch("TELESCOPE_PARK")
        while(not(unpark)):
            unpark=self.getSwitch("TELESCOPE_PARK")
        unpark[0].setState(PyIndi.ISS_OFF)
        unpark[1].setState(PyIndi.ISS_ON)
        self.indiclient.sendNewSwitch(unpark)