import PyIndi
import indi
class BaseDevice:
    device: PyIndi.BaseDevice
    name = ""
    indiclient: indi.IndiClient
    def __init__(self, name, device, indiclient) -> None:
        self.name = name
        self.device = device
        self.indiclient = indiclient

    def Connect(self):
        print("connecting to device "+ self.name)
        if (not(self.device.isConnected())):
            print("device is not connected")
            device_connect = self.device.getSwitch("CONNECTION")
            device_connect[0].s = PyIndi.ISS_ON
            device_connect[1].s = PyIndi.ISS_OFF
            self.indiclient.sendNewSwitch(device_connect)


    def Disconnect(self):
        if (self.device.isConnected()):
            print("device is connected")
            device_connect = self.device.getSwitch("CONNECTION")
            device_connect[0].s = PyIndi.ISS_OFF
            device_connect[1].s = PyIndi.ISS_ON
            self.indiclient.sendNewSwitch(device_connect)

    def getSwitch(self, switchName):
        return self.device.getSwitch(switchName)
            
    def getNumber(self, numberName):
        return self.device.getNumber(numberName)
    
    def getText(self, text):
        return self.device.getText(text)
    
    def sendNewSwitch(self, switch):
        self.indiclient.sendNewSwitch(switch)  

    def sendNewNumber(self, number):
        self.indiclient.sendNewNumber(number)        
    
    def sendNewText(self, text):
        self.indiclient.sendNewText(text)