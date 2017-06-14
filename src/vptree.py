import numpy as np

class vptree:
    def __init__(self, data):
        self.head = None
        self.size = data.shape[1]
        self.data = None

    def sort(self):
        self.data = np.sort(self.data, 
    def construct(self, data):
        self.data = data

        data_az = data[0]
        data_po = data[1]

        head_az = 

        
    def metric(self, az_a, po_a, az_b, po_b):
        return np.arccos(np.sin(po_a)*np.sin(po_b)*np.cos(az_a-az_b) + np.cos(po_a)*np.cos(po_b)

    
