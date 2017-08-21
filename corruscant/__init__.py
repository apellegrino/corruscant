from os.path import abspath, dirname
from ctypes import CDLL

import os
import re

PROJECT_PATH = dirname(abspath(__file__))

def loadlib(name):
    pattern = re.compile(name)
    binpath = os.path.join(PROJECT_PATH, "bin")
    for filename in os.listdir(binpath):
        m = re.match(pattern, filename)
        if m:
            return CDLL(os.path.join(binpath, filename))
    else:
        raise OSError("Shared library not found. Did you remember to \"make\"?")

kdlib = loadlib("libkdtree")
coordlib = loadlib("libcoords")
