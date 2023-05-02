import PyIndi
import time
import sys
import device.mount
import device.ccd
import device.base_device
import polar_align.polar_align_service
import indi
import traceback

# connect the server
indiclient=indi.IndiClient()
indiclient.setServer("localhost",7624)

if (not(indiclient.connectServer())):
        print("No indiserver running on "+indiclient.getHost()+":"+str(indiclient.getPort())+" - Try to run")
        print("  indiserver indi_simulator_telescope indi_simulator_ccd")
        sys.exit(1)

# connect the scope
#telescope="Telescope Simulator"
telescope="EQMod Mount"
device_telescope=None
telescope_connect=None

# get the telescope device
device_telescope=indiclient.getDevice(telescope)
max=0
while not(device_telescope) or max < 3:
    print(f"connecting telescope, try {max}")
    time.sleep(3)
    device_telescope=indiclient.getDevice(telescope)
    max+=1
print(device_telescope)
if not(device_telescope):
    print("no telescope found, exiting...")
    sys.exit(1)
mount=device.mount.Mount(telescope, device_telescope, indiclient)
mount.Connect()
telescope_info=mount.getNumber("TELESCOPE_INFO")
telescope_info[0].setValue(73)
telescope_info[1].setValue(400)
mount.sendNewNumber(telescope_info)
mount.unpark()
ccd="ZWO Camera"
ccd="CCD Simulator"
device_ccd=indiclient.getDevice(ccd)
while(not(device_ccd)):
    print("retrying")
    time.sleep(3)
    device_ccd=indiclient.getDevice(ccd)
ccdDevice=device.ccd.Ccd(ccd, device_ccd, indiclient)

while(not(ccdDevice.device.isConnected())):
    print("retrying")
    ccdDevice.Connect()
    time.sleep(3)

indiclient.setBLOBMode(PyIndi.B_ALSO, ccdDevice.name, "CCD1")
uploadDir=ccdDevice.getText("UPLOAD_SETTINGS")
while(not(uploadDir)):
    uploadDir=ccdDevice.getText("UPLOAD_SETTINGS")
uploadDir[0].text="/home/nicolas/indiPhotos"
uploadDir[1].text="test_"
print("uploadDir ", uploadDir[0].text)
indiclient.sendNewText(uploadDir)

uploadMode=ccdDevice.getText("UPLOAD_MODE")
while(not(uploadMode)):
    uploadMode=ccdDevice.getSwitch("UPLOAD_MODE")
uploadMode[0].setState(PyIndi.ISS_OFF)
uploadMode[1].setState(PyIndi.ISS_OFF)
uploadMode[2].setState(PyIndi.ISS_ON)
indiclient.sendNewSwitch(uploadMode)
indiclient.addBlobHandler(ccdDevice)
ccdDevice.setActiveDevices(device=device_telescope)
polar_align = polar_align.polar_align_service.PolarAlignService(mount, ccdDevice, indiclient)
polar_align.StartPolarAlign(None, None)
try:
    for i in range(10):
        polar_align.ComputePolarError()
except:
    print(traceback.print_exc())
mount.Park()