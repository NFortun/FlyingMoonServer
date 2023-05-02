import PyIndi
import threading

class BlobNotifier:
    blobEvent: threading.Event
    def notifyBlob(self):
        print("notify blob")
        self.blobEvent.set()

class IndiClient(PyIndi.BaseClient):
    blobHandlers = []
    def __init__(self):
        super(IndiClient, self).__init__()
    def newDevice(self, d):
        pass
    def newProperty(self, p):
        pass
    def removeProperty(self, p):
        pass
    def newBLOB(self, bp):
        print("new BLOB ", bp.name)
        for handler in self.blobHandlers:
            handler.notifyBlob()
    def newSwitch(self, svp):
        pass
    def newNumber(self, nvp):
        pass
    def newText(self, tvp):
        pass
    def newLight(self, lvp):
        pass
    def newMessage(self, d, m):
        pass
    def serverConnected(self):
        pass
    def serverDisconnected(self, code):
        pass
    def addBlobHandler(self, handler):
        self.blobHandlers.append(handler)
