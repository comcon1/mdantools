#!/usr/bin/python
import re

class GmxIndex:
  def __init__(self, filename):
    f = open(filename, 'r')
    self._data = {}
    nm = re.compile('\[\s*([a-z\+0-9\-_]+)\s*\]', re.IGNORECASE) #"[ SOL ]"
    # nums = re.compile( TODO
    curnm = ''
    for line in f:
      m = nm.search(line)
      if m != None:
        curnm = m.group(1)
        self._data[curnm] = []
        continue
      else:
        ar = map(int,line.split())
        self._data[curnm] += ar
    f.close()

  def getNdx(self, name):
    return self._data[name]

  def getRNdx(self, name):
    return map(lambda x: x-1, self._data[name])
