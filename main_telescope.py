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
while not(device_telescope) and max < 3:
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

#for p in device_telescope.getProperties():
#    print(p.getName())
#telescope_info=mount.getNumber("TELESCOPE_INFO")
#telescope_info[0].setValue(73)
#telescope_info[1].setValue(400)
#mount.sendNewNumber(telescope_info)
mount.unpark()
for p in device_telescope.getProperties():
    if (p.getName() == "TELESCOPE_PARK"):
         park=device_telescope.getSwitch("TELESCOPE_PARK")
         print(park[0].getName()+" "+ str(park[0].getState()))
         print(park[1].getName()+" "+ str(park[1].getState()))
         
         
mount.goto(19,45)
time.sleep(5)
mount.goto(15)
time.sleep(5)
mount.Park()
#for p in device_telescope.getProperties():
#    if (p.getName() == "TELESCOPE_PARK"):
#         park=device_telescope.getSwitch("TELESCOPE_PARK")
#         print(park[0].getName()+" "+ str(park[0].getState()))
#         print(park[1].getName()+" "+ str(park[1].getState()))
#         
   