# Object used to duplicate sys.stdout and sys.stderr to a log file. Emulates
# tee behaviour. Copied from
# http://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-stdout-to-a-log-file-in-python

import sys


class Tee(object):

    def __init__(self, name, mode):
        self.file = open(name, mode, 0)  # '0' is for no buffering
        self.fileno = self.file.fileno
        self.flush = self.file.flush
        self.stdout = sys.stdout
        self.stderr = sys.stderr
        sys.stdout = self
        sys.stderr = self

    def __del__(self):
        sys.stdout = self.stdout
        sys.stderr = self.stderr
        self.file.close()

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)
