
import yoda, rivet, os
from math import sqrt

inFile  = 'Rivet.yoda'

hists = yoda.read( inFile )
tags = sorted(hists.keys())

# probably more fancy than it needs to be ...
def getRivetRefData(anas=None):
    "Find all Rivet reference data files"
    refhistos = {}
    for var in ("RIVET_ANALYSIS_PATH", "RIVET_DATA_PATH", "RIVET_REF_PATH", "RIVET_INFO_PATH", "RIVET_PLOT_PATH"):
        if var in os.environ:
            abspaths = map(os.path.abspath, os.environ[var].split(":"))
            os.environ[var] = ":".join(abspaths)
    rivet_data_dirs = rivet.getAnalysisRefPaths()
    dirlist = [ ]
    for d in rivet_data_dirs:
        import glob
        if anas is None:
          dirlist.append(glob.glob(os.path.join(d, '*.yoda')))
        else:
          for a in anas:
            dirlist.append(glob.glob(os.path.join(d, a+'*.yoda')))
    for filelist in dirlist:
        for infile in filelist:
            analysisobjects = yoda.read(infile)
            for path, ao in analysisobjects.iteritems():
                aop = rivet.AOPath(ao.path)
                if aop.isref():
                    ao.path = aop.basepath(keepref=False)
                    refhistos[ao.path] = ao
    return refhistos

# get hold of relevant objects in reference data files
refhistos = getRivetRefData(['ATLAS_2017_I1624693'])

def constructDiff(hist):
    '''This function produces a (data - MC)/sigma version of the Dalitz (2D) plot.'''
    path_data = hist.annotation('Path') # data central value
    path_stat = hist.annotation('Path').replace('d03', 'd04') # data statistical uncertainty
    path_unco = hist.annotation('Path').replace('d03', 'd05') # data uncorrelated uncertainty
    data = refhistos[path_data]
    stat = refhistos[path_stat]
    unco = refhistos[path_unco]
    data_integral = sum([ p.z for p in data.points ])
    hist.normalize(data_integral) 
    rtn = hist.clone();  rtn.reset()
    mc = yoda.mkScatter(hist)
    for i in range(rtn.numBins):
        sigma = sqrt(stat.points[i].z ** 2 + unco.points[i].z ** 2 + (mc.points[i].zErrs[0]) ** 2)
        newz = (data.points[i].z - mc.points[i].z) / sigma if sigma else 0.0
        #newz = (0.01 * mc.points[i].z - data.points[i].z) / sigma if sigma else 0.0
        rtn.fillBin(i, newz)
    return rtn

# this is where the magic happens
f = open('%s_processed.yoda' % inFile[:-5], 'w')
for h in tags:
    if 'd03-x01-y01' in h:
        hdiff = constructDiff(hists[h])
        outName = h.replace('y01', 'y02')
        hdiff.setAnnotation('Path', outName)
        yoda.writeYODA(hdiff, f)
    yoda.writeYODA(hists[h], f)
f.close()

