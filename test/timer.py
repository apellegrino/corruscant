from time import clock

class timer:
    def __init__(self):
        self.time = 0.0
        self.started = 0.0

    def start(self):
        self.started = clock()

    def stop(self):
        self.time += (clock() - self.started)

    def read(self):
        print "Elapsed Time: %f s" % self.time
        
