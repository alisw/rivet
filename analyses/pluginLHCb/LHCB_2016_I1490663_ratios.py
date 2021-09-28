#! /usr/bin/env python

# YODA script to generate the ratio plots between 13 (LHCB_2015_I1396331) and 5 TeV (LHCB_2016_I1490663) c.o.m. energies
# for D0, D+, Ds+ and D*(2010)+ and their charge conjugates (in this order)

import yoda
import os
import re
import sys


aname05 = 'LHCB_2016_I1490663'
aname13 = 'LHCB_2015_I1396331'
nybins = 5


def showHelp():
    print('usage:')
    print('    <python-script> <histo-numer.yoda> <histo-denom.yoda>')


def makeRatio(yd, sref, hn, hd):
    nedges = list(sref.xMins())
    nedges.append(sref.xMaxs()[-1])
    # re-bin to reference data edges
    hn.rebinTo(nedges)
    hd.rebinTo(nedges)
    tt = os.path.split(hd.path())[0]
    tt = tt.replace('/RAW/', '/')   # put Scatter2D object in root
    tpath = unicode(os.path.join(tt, os.path.split(sref.path())[1]))
    yd[tpath] = yoda.core.divide(hn, hd)
    # ensure Scatter2D object saved to correct path (!)
    yd[tpath].setPath(tpath)


# start code
if __name__ != "__main__":
    raise NotImplementedError("Please, run this file as a script. It is not a module to import.")

if len(sys.argv) < 3:
    showHelp()
    sys.exit(1)

if any([not (x.endswith('.yoda') and os.path.isfile(x)) for x in sys.argv[1:2]]):
    showHelp()
    sys.exit(1)

arg13 = sys.argv[1]
arg05 = sys.argv[2]

# first entry in tuple is ratio dataset id, second is denominator dataset id and third is the numerator dataset id
ratiodef = [('d05', 'd01', 'd01'), ('d06', 'd02', 'd02'), ('d07', 'd03', 'd03'), ('d08', 'd04', 'd04')]
pre13 = [re.compile(r'/RAW/%s/%s-x01-y0[1-5]' % (aname13, ds[1])) for ds in ratiodef]
pre05 = [re.compile(r'/RAW/%s/%s-x01-y0[1-5]' % (aname05, ds[2])) for ds in ratiodef]
ref05 = [re.compile(r'/REF/%s/%s-x01-y0[1-5]' % (aname05, ds[0])) for ds in ratiodef]
d13 = yoda.read(arg13, patterns=pre13)
if len(d13.values()) != nybins*len(ratiodef) or aname13 not in d13.values()[0].path():
    print('ERR: Wrong or incomplete 13 TeV data set.')
    showHelp()
    sys.exit(2)
d05 = yoda.read(arg05, patterns=pre05)
if len(d05.values()) != nybins*len(ratiodef) or aname05 not in d05.values()[0].path():
    print('ERR: Wrong or incomplete 5 TeV data set.')
    showHelp()
    sys.exit(2)

ys = yoda.read(arg05)
srefs = yoda.read('%s.yoda' % aname05, patterns=ref05)
for i in range(1, nybins+1):
    xyset = '-x01-y0%d' % i
    for rdef in ratiodef:
        rk05 = u'/REF/%s/%s%s' % (aname05, rdef[0], xyset)
        dk13 = u'/RAW/%s/%s%s' % (aname13, rdef[1], xyset)
        dk05 = u'/RAW/%s/%s%s' % (aname05, rdef[2], xyset)
        makeRatio(ys, srefs[rk05], d13[dk13], d05[dk05])

# use denominator histo file name when saving, yet avoiding overwrite
bname = os.path.splitext(arg05)[0]
rFile = '%s+ratios.yoda' % bname
yoda.writeYODA(ys,  rFile)

sys.exit(0)
