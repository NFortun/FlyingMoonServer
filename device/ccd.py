import os
import time
import PyIndi
import threading
from device.base_device import BaseDevice
from indi import BlobNotifier


class Ccd(BaseDevice, BlobNotifier):
    def __init__(self, name, mount_device, indiclient) -> None:
        super().__init__(name, mount_device, indiclient)

    def setExposure(self, exposure):
        ccd_exposure = self.getNumber("CCD_EXPOSURE")
        ccd_exposure[0].setValue(exposure)
        self.sendNewNumber(ccd_exposure)
    
    def setActiveDevices(self, device: PyIndi.BaseDevice):
        active_devices = self.getText("ACTIVE_DEVICES")
        active_devices[0].setText(device.getDeviceName())
        self.sendNewText(active_devices)

    def takeExposure(self) -> str: 
        print("taking exposure")
        uploadDir=self.getText("UPLOAD_SETTINGS")
        now = time.time()
        uploadDir[1].text=str(now)
        self.indiclient.sendNewText(uploadDir)
        ccd_ccd1=self.device.getBLOB("CCD1")
        self.blobEvent=threading.Event()
        self.blobEvent.wait()
        for blob in ccd_ccd1:
            print("name: ", blob.name," size: ", blob.size," format: ", blob.format)
            # pyindi-client adds a getblobdata() method to IBLOB item
            # for accessing the contents of the blob, which is a bytearray in Python
            fits=blob.getblobdata()
            print("fits data type: ", type(fits))
        return os.path.join(uploadDir[0].text, uploadDir[1].text)