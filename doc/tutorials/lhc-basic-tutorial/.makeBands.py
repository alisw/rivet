#!/usr/bin/env python

import yoda
import numpy as np

FILE_NAME = 'Rivet.yoda'

SCALES = [
 'MUR0.5_MUF0.5_PDF261000_PSMUR0.5_PSMUF0.5',
 'MUR0.5_MUF1_PDF261000_PSMUR0.5_PSMUF1',
 'MUR1_MUF0.5_PDF261000_PSMUR1_PSMUF0.5',
 'MUR1_MUF1_PDF261000',
 'MUR1_MUF2_PDF261000_PSMUR1_PSMUF2',
 'MUR2_MUF1_PDF261000_PSMUR2_PSMUF1',
 'MUR2_MUF2_PDF261000_PSMUR2_PSMUF2',
]

def constructBand():
  # remove .yoda from file name
  fName = FILE_NAME[:-5] if FILE_NAME.endswith('.yoda') else FILE_NAME
  # open output files
  outErrStat = open('%s_stat_only.yoda' % fName, 'w')
  outErrStatME = open('%s_stat_MEonly.yoda' % fName, 'w')
  outErrStatMEPS = open('%s_full_band.yoda' % fName, 'w')
  # read in input files
  hists = yoda.read('%s.yoda' % fName)
  # loop over histograms
  for name in hists:
    # skip auxiliary objects
    if 'RAW' in name or name.startswith('/_'):
      continue
    # only worry about Histo1D objects for now:
    if 'Histo1D' not in str(hists[name]):
      continue
    # skip variation objects
    if name.endswith(']'):
      continue
    # now only process nominal objects, check:
    print ('Processing: %s' % name)
    # stats-only is trivial -- copy goes
    # straight into respective output file:
    yoda.writeYODA(hists[name].mkScatter(), outErrStat) 
    binw  = np.array([ b.xWidth() for b in hists[name].bins() ]) # bin widths
    noms  = np.array([ b.sumW()   for b in hists[name].bins() ]) # nominal central values
    stat2 = np.array([ b.sumW2()  for b in hists[name].bins() ]) # squared statistical uncertainties
    # calculate maximum-shift envelope from 7-point scale variations
    scaleMEPSup = np.array(noms);  scaleMEup = np.array(noms)
    scaleMEPSdn = np.array(noms);  scaleMEdn = np.array(noms)
    #-----------------------------------------------------------------------------------
    # TODO: Calculate the envelope of the 7-point scale variations
    # HINT: If you loop over the strings in SCALES (defined near the top of this file), 
    #       you can use hists["%s[%s]" % (name, string_in_SCALES)] to retrieve the
    #       respective ME+PS scale variation histogram from the dictionary "hists"
    #       and hists["%s["ME_ONLY_"%s]" % (name, string_in_SCALES)] to retrieve the
    #       corresponding ME-only version.
    # Add your code here:
    for scale in SCALES:
      var = np.array([ b.sumW() for b in hists["%s[%s]" % (name, scale)].bins() ])
      scaleMEPSup = list(map(max, zip(scaleMEPSup, var)))
      scaleMEPSdn = list(map(min, zip(scaleMEPSdn, var)))
      var = np.array([ b.sumW() for b in hists["%s[ME_ONLY_%s]" % (name, scale)].bins() ])
      scaleMEup   = list(map(max, zip(scaleMEup, var)))
      scaleMEdn   = list(map(min, zip(scaleMEdn, var)))
    #-----------------------------------------------------------------------------------
    # TODO: Add the respective scale-uncertainty component in quadrature with 
    #       the statistical uncertainty. Don't forget to divide by bin width.
    # Add your code here:
    systMEPSup = np.sqrt(stat2 + (scaleMEPSup - noms) ** 2) / binw
    systMEPSdn = np.sqrt(stat2 + (noms - scaleMEPSdn) ** 2) / binw
    systMEup = np.sqrt(stat2 + (scaleMEup - noms) ** 2) / binw
    systMEdn = np.sqrt(stat2 + (noms - scaleMEdn) ** 2) / binw
    #-----------------------------------------------------------------------------------
    hMEPS = hists[name].mkScatter()
    hME = hMEPS.clone()
    for i in range(hMEPS.numPoints()):
      hMEPS.point(i).setYErrs(systMEPSdn[i], systMEPSup[i])
      hME.point(i).setYErrs(systMEdn[i], systMEup[i])
    # save in respective output file and continue
    yoda.writeYODA(hME, outErrStatME)
    yoda.writeYODA(hMEPS, outErrStatMEPS)
  outErrStat.close()
  outErrStatME.close()
  outErrStatMEPS.close()


if __name__ == '__main__':
  constructBand()
