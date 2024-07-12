from __future__ import print_function
import rivet, yoda
import os, glob, logging, re
from math import sqrt
from rivet.plotting import plot2yaml
from rivet.plotting.conversion_tools import type_conversion


# TODO: add more descriptive docstrings to all functions.

def _sanitise_string(s):
    s = s.replace('#','\\#')
    s = s.replace('%','\\%')
    return s


def _parse_args(args):
    """Look at the argument list and split it at colons, in order to separate
    the file names from the plotting options. Store the file names and
    file specific plotting options.

    Parameters
    ----------
    args : list[str]
        List of arguments which were previously passed to `rivet-cmphistos`.
        Format will be ['filename.yoda:key=value', ..., 'PLOT:key=value:key=value']

    Returns
    -------
    filelist : list[str]
        Raw names of the files, i.e. the first part of each string in args.
    filenames : list[str]
        Names of the files. If Name=value is passed as a plot option after a file name, this will become the filename. Otherwise, it will use the same value as in filelist.
    plotoptions : dict[str, dict[str, str]]
        Dictionary of plot options.
        The key will be the file name (i.e., same value as in filenames) and the value will be a dict of strings with plot options.
        One of the keys will also be PLOT (if it was passed in as an argument to args, which contains all plot options that will be applied to the entire figure.

    Examples
    --------
    >>> _parse_args(['mc1.yoda:Title=example title:Name=example name 1', 'mc2.yoda', 'PLOT:LogX=1'])
    (['mc1.yoda', 'mc2.yoda'],
     ['example name 1', 'mc2.yoda'],
     {
         'example name 1': {'Title': 'example title', 'Name': 'example name 1'},
         'mc2.yoda': {'Title': 'mc2'},
         'PLOT': {'LogX': '1'}
    })

    Note
    ----
    Some matplotlib line styles contain ':', which would not work with current code. TODO: change delimiter?
    """
    # TODO: remove filenames since they exist as keys in plotoptions?
    filelist = []
    filenames = []
    plotoptions = {}
    for a in args:
        asplit = a.split(':')
        path = asplit[0]
        if path != "PLOT":
            filelist.append(path)
            filenames.append(path)
        plotoptions[path] = {}
        has_title = False
        has_name = ""
        for i in range(1, len(asplit)):
            ## Add 'Title' if there is no = sign before math mode
            if '=' not in asplit[i] or ('$' in asplit[i] and asplit[i].index('$') < asplit[i].index('=')):
                asplit[i] = 'Title=%s' % asplit[i]
            if asplit[i].startswith('Title='):
                has_title = True
            key, value = asplit[i].split('=', 1)
            plotoptions[path][key] = type_conversion(value)
            if asplit[i].startswith('Name=') and path != "PLOT":
                has_name = asplit[i].split('=', 1)[1]
                filenames[-1] = has_name
        if has_name != "":
            plotoptions[has_name] = plotoptions[path]
            del plotoptions[path]
        if path != "PLOT" and not has_title:
            plotoptions[has_name if has_name != "" else path]['Title'] = _sanitise_string(os.path.basename( os.path.splitext(path)[0] ))
    return filelist, filenames, plotoptions


def _get_histos(filelist, filenames, plotoptions, path_patterns = [], path_unpatterns = []):
    """Loop over all input files. Only use the first occurrence of any REF-histogram
    and the first occurrence in each MC file for every MC-histogram."""

    refhistos, mchistos = {}, {}
    for infile, inname in zip(filelist, filenames):
        mchistos.setdefault(inname, {})
        try:
            analysisobjects = yoda.read(infile, patterns=path_patterns, unpatterns=path_unpatterns)
        except IOError as e:
            print("File reading error:", e.strerror)
            sys.exit(1)
        for path, ao in analysisobjects.items():

            # Make a path object and ensure the path is in standard form.
            try:
                aop = rivet.AOPath(path)
            except Exception as e:
                print("Found analysis object with non-standard path structure:", path, "... skipping")
                continue

            ## We don't plot data objects with path components hidden by an underscore prefix
            if aop.istmp() or aop.israw():
                continue

            # Convert non-scatter objects to scatter
            if "Scatter" not in ao.type():
                ao = ao.mkScatter()

            ## Add it to the ref or mc paths, if this path isn't already known
            basepath = aop.basepath(keepref=False)
            defaultWeightName = plotoptions[inname].get('DefaultWeight', '0')
            if aop.isref() and basepath not in refhistos:
                ao.setPath(aop.varpath(keepref=False, defaultvarid=defaultWeightName))
                refhistos[basepath] = ao
            else: #if basepath not in mchistos[infile]:
                mchistos[inname].setdefault(basepath, {})[aop.varid(defaultWeightName)] = ao

    return refhistos, mchistos


def _get_rivet_ref_data(anas, path_patterns, path_unpatterns):
    """Find all Rivet reference data files"""
    refhistos = {}
    rivet_data_dirs = rivet.getAnalysisRefPaths()
    dirlist = []
    for d in rivet_data_dirs:
        if anas is None:
            dirlist.append(glob.glob(os.path.join(d, '*.yoda*')))
        else:
            #dirlist.append([os.path.join(d, a+'.yoda*') for a in anas])
            for a in anas:
                res = glob.glob(os.path.join(d, a+'.yoda'))
                if len(res) == 0:
                    res = glob.glob(os.path.join(d, a+'.yoda.gz'))
                if len(res) != 0:
                    dirlist.append(res)
    for filelist in dirlist:
        # TODO: delegate to _get_histos?
        for infile in filelist:
            analysisobjects = yoda.read(infile, patterns=path_patterns, unpatterns=path_unpatterns)
            for path, ao in analysisobjects.items():
                aop = rivet.AOPath(ao.path())
                if aop.isref():
                    ao.setPath(aop.basepath(keepref=False))
                    refhistos[ao.path()] = ao
    return refhistos


def get_nominal_key(listOfHistoKeys):
    """try to find the key corresponding to the nominal histogram,
       which seems to differ between different YODA files
    """
    name = '0'
    if 'nominal' in listOfHistoKeys: name='nominal'
    elif 'yoda' in listOfHistoKeys:  name='yoda'
    elif '0' in listOfHistoKeys:  name='0'
    return name


def _make_output(plot_id, plotdirs, config_files, mchistos, refhistos, reftitle, plotoptions,
                 style, rc_params, mc_errs, nRatioTicks, skipWeights, removeOptions, deviation, canvasText, verbose = False):
    """Create output dictionary for the plot_id.

    Parameters
    ----------
    plot_id : str
        ID, usually of the format AnalysisID/HistogramID.
    plotdirs : list[str]
        All directories to look for .plot files at.
    config_files : list[str]
        Additional plot settings that will be applied to all figures.
    mchistos : dict
        Dictionary of the Monte Carlo YODA histograms.
        The structure is {filename: {plot_id: {"0": yoda_histogram1, "1": yoda_histogram2, ...}}}
        Usually only "0" exists as the innermost key.
    refhsitos : dict
        Dictionary of the reference analysis data YODA histograms.
    plotoptions : dict[str, dict[str, str]]
        Dict containing all plot options for all histograms and all plots.
    mc_errs : bool
        See rivet_mkdat
    style : str
        A predefined name of a style.
    removeOptions : bool
        If true, prevents appending the options string to the legend label
    deviation : bool
        If true, express compatability between curve and ref. data in terms 
        of standard deviations in ratio panel.
    rc_params : dict[str, str]
        Dict of rcParams that will be added to the rcParams section of the output .dat file.

    Returns
    -------
    outputdict : dict
        Correctly formatted dictionary that can be passed to `yaml.dump` to write to an output file.
    """
    outputdict = {}
    plot_configs = plot2yaml.get_plot_configs(plot_id, plotdirs=plotdirs, config_files=config_files)
    outputdict['plot features'] = plot_configs
    outputdict['plot features']['Deviation'] = deviation
    
    # only write extra info to the .dat file if specified by user
    if nRatioTicks !=1: outputdict['plot features'].update({"nRatioTicks": nRatioTicks})
    if canvasText != None: outputdict['plot features'].update({"canvasText" : canvasText})
    outputdict['plot features'].update(plotoptions.get('PLOT', {}))
    outputdict['rcParams'] = rc_params
    outputdict['style'] = style
    outputdict['stylepath'] = '../'
    outputdict['histograms'] = {}

    componentNames = ['BandComponentPDF', 'BandComponentEnv']

    if plot_id in refhistos and refhistos[plot_id].type() != 'Scatter1D':
        refhistos[plot_id].setAnnotation('IsRef', True)
        outputdict['histograms'][reftitle] = {'nominal': refhistos[plot_id]} # this is where ErrorBreakdown is included?
        outputdict['histograms'][reftitle]['IsRef'] = True
        outputdict['histograms'][reftitle]['Title'] = 'Data'
        outputdict['plot features']['RatioPlotYLabel'] = 'MC/Data'
        outputdict['plot features']['RatioPlot'] = True
    lhapdfCheck = True
    for filename, mchistos_in_file in mchistos.items():
        for plot_id_with_anaopt in sorted(mchistos_in_file):
            if rivet.stripOptions(plot_id_with_anaopt) != plot_id:
                continue
            histogroup = mchistos_in_file[plot_id_with_anaopt]
            # TODO Counters/Scatter1D currently not supported
            if next(iter(histogroup.values())).type() == 'Scatter1D':
                continue

            label = rivet.extractOptionString(plot_id_with_anaopt)
            outputdict['histograms'][filename+label] = {}

            thisFilePlotOptions = dict(plotoptions.get(filename, {}))
            # add options string to legend entry
            newtitle = thisFilePlotOptions.get('Title', '')
            if not removeOptions:
                newtitle += label
            thisFilePlotOptions['Title'] = newtitle
            outputdict['histograms'][filename+label].update(thisFilePlotOptions)

            makePDFBand  = thisFilePlotOptions['BandComponentPDF'] \
                           if 'BandComponentPDF' in thisFilePlotOptions else ''
            makeEnvelope = thisFilePlotOptions['BandComponentEnv'] \
                           if 'BandComponentEnv' in thisFilePlotOptions else ''

            # check if lhapdf is available
            if makePDFBand and lhapdfCheck:
                try:
                    import lhapdf
                    lhapdf.setVerbosity(0)
                    lhapdfCheck = False
                except ImportError as e:
                    print("LHAPDF not available! Need this to construct PDF band:",
                            f" failing `import {e.name}`")
                    exit(1)

            nominalVariationKey = get_nominal_key(mchistos_in_file[plot_id_with_anaopt].keys())
            if nominalVariationKey == None:
                raise NameError("Could not find nominal variation weight!")
            obj_type = mchistos_in_file[plot_id_with_anaopt][nominalVariationKey].type()

            nomVals = None
            pdf_matches = { }; env_matches = { }
            PDFvars = [ [] for _ in makePDFBand.split() ]
            PDFsets = [ None for _ in makePDFBand.split() ]
            Enverrors = [ [] for _ in makeEnvelope.split() ]
            for histogramkey, histogram in histogroup.items():
                isNominal = (nominalVariationKey == histogramkey)
                # Maybe add this mc_errs option to the plotoptions dict and only
                # pass the plotoptions dict to the function?
                outputdict['histograms'][filename+label]['ErrorBars'] = mc_errs

                thisObj = histogram.mkScatter() if 'Histo1D' in obj_type else histogram
                # central values of current object
                central_values = [ p.y() for p in thisObj.points() ]

                if isNominal:
                    nominalScatter = histogram
                    outputdict['histograms'][filename+label]['nominal'] = histogram
                    nomvals = list(central_values)

                for i, prescription in enumerate(makePDFBand.split()):
                    if verbose and prescription not in pdf_matches:
                      pdf_matches[prescription] = [ ]
                    if any([ re.search(pat, histogramkey) for pat in prescription.split(',') ]):
                        if verbose:
                            pdf_matches[prescription].append(histogramkey)

                        # store values from PDF variation
                        PDFvars[i].append(central_values)

                        # initialize pdf set object from lhapdf
                        if PDFsets[i] is None:
                            lhapdfID = int(re.search('PDF[0-9]*', histogramkey).group(0)[3:])
                            PDFsets[i] = lhapdf.mkPDF(lhapdfID).set()

                for i, prescription in enumerate(makeEnvelope.split()):
                    if verbose and prescription not in env_matches:
                      env_matches[prescription] = [ ]
                    if isNominal or any([ re.search(pat, histogramkey) for pat in prescription.split(',') ]):
                        if verbose:
                            env_matches[prescription].append(histogramkey)
                        if not Enverrors[i]:
                            Enverrors[i] = [ list(central_values), list(central_values) ]
                        else:
                            Enverrors[i][0] = list(map(min, zip(Enverrors[i][0], central_values)))
                            Enverrors[i][1] = list(map(max, zip(Enverrors[i][1], central_values)))

                # don't plot multiweights if already plotting a band
                if not skipWeights and not isNominal and not makeEnvelope and not makePDFBand:
                    if not rivet.extractWeightName(plot_id_with_anaopt).startswith('EXTRA'):
                        outputdict['histograms'][filename+label]['multiweight'+histogramkey] = histogram

            if verbose:
                for pat in pdf_matches:
                  print ("PDF prescription \"%s\" matches:" % pat)
                  print (pdf_matches[pat])
                for pat in env_matches:
                  print ("Envelope prescription \"%s\" matches:" % pat)
                  print (env_matches[pat])
                del pdf_matches, env_matches

            PDFerrors = [ ]
            for pdf_set, pdf_vars in zip(PDFsets, PDFvars):
                # if number of PDF variations if off by 1,
                # probably needs the nominal
                if len(pdf_vars) == int(pdf_set.size) - 1:
                    pdf_vars.append(nomVals)
                elif len(pdf_vars) != int(pdf_set.size):
                    raise ValueError("Number of matched PDF variations is %i, expected %s!" % (len(pdf_vars), pdf_set.size))
                pdf_vars = list(map(list,zip(*pdf_vars))) # transpose
                try:
                    # calculate uncertainties from all PDFs multiweight
                    # histos that matched regex from the user
                    uncertainties = [ pdf_set.uncertainty(binVars) for binVars in pdf_vars ]
                    PDFerrors.append([ (unc.errminus, unc.errplus) for unc in uncertainties ])

                except RuntimeError:
                    print("Error in constructing the PDFset. Skipping.")

            # let user ask for a band, if no BandComponentEnv/PDF provided
            # this will just be a band with stat. errors
            if makePDFBand or makeEnvelope:

                BandScatter = nominalScatter.clone()
                # iterate over bins
                for ibin, y in enumerate(nominalScatter.yVals()):
                    totErrDn, totErrUp = BandScatter.point(ibin).yErrs()
                    totErrDn = totErrDn*totErrDn
                    totErrUp = totErrUp*totErrUp

                    # add PDF uncertainty in quadrature
                    for errs in PDFerrors:
                        totErrDn += errs[ibin][0]*errs[ibin][0]
                        totErrUp += errs[ibin][1]*errs[ibin][1]

                    # add Envelope uncertainty in quadrature
                    for errDn, errUp in Enverrors:
                        absEnvDn = y - errDn[ibin]
                        absEnvUp = errUp[ibin] - y
                        totErrDn += absEnvDn*absEnvDn
                        totErrUp += absEnvUp*absEnvUp

                    # Scatter object with total Band uncertainty
                    BandScatter.point(ibin).setYErrs(sqrt(totErrDn), sqrt(totErrUp))
                outputdict['histograms'][filename+label]['BandUncertainty'] = BandScatter

    # Remove all sections of the output_dict that do not contain any information.
    # A list of keys is first created. Otherwise, it will raise an error since the size of the dict changes.
    dict_keys = list(outputdict.keys())
    for key in dict_keys:
        if not outputdict[key]:
            del outputdict[key]
    return outputdict


def assemble_plotting_data(args, path_pwd=True, reftitle='Data', rivetrefs=True,
                           path_patterns=[], path_unpatterns=[],
                           plotinfodirs=[], style='default', config_files=[],
                           hier_output=False, outdir='.', mc_errs=True,
                           rivetplotpaths=True, analysispaths=[], verbose=False,
                           writefiles=False, nRatioTicks=1, skipWeights=False, 
                           removeOptions = False, deviation=False,                    
                           canvasText=None):
    """Create a dictionary of the plotting data that can be turned
    into self-consistent Python executables.

    Parameters
    ----------
    args : Iterable[str]
        Non-keyword arguments that were previously passed to rivet-cmphistos.
        E.g., ['mc1.yoda', 'mc2.yoda:Title=example title', 'PLOT:LogX=1']
    path_pwd : bool
        Search for plot files and reference data files in current directory.
    reftitle : str
        Legend name of the reference data in the plots.
    rivetrefs : bool
        If False, don't use Rivet reference data files
    path_patterns : Iterable[str]
        Only write out histograms whose $path/$name string matches these regexes.
        The argument may also be a text file.
    path_unpatterns : Iterable[str]
        Exclude histograms whose $path/$name string matches these regexes
    plotinfodirs : list[str]
        Directory which may contain plot header information (in addition to standard Rivet search paths).
    style : str
        Set the style of all plots and additional rcParams.
        Format is style:key=value:key2=value2...
        The first part of the string must be a name of a builtin style (e.g. 'default').
        The other keys and values must be valid rcParams.
        However, the validity is not checked by this function.
    config_files : list[str]
        Additional plot config file(s).
        Settings will be included in the output configuration.
        ~/.make-plots will automatically be added.
    hier_output : bool
        Write output .dat files into a directory hierarchy which matches the analysis paths.
    outdir : str
        Write dat files into this directory.
    mc_errs : bool
        If True, add the errors of the Monte-Carlo histograms.
    rivetplotpaths : bool
        Search for .plot files in the standard Rivet plot paths.
    verbose : bool
        If True, write more information to stdout.
    writefiles : bool
        If True, write the created dicts to dat files.
        This is used if one wants the intermediate format for later use or if one only calls this function and not rivet-mkhtml.
    nRatioTicks: int
        Number of minor ticks between major ticks, can be specified in rivet-mkhtml
    deviation: bool
        Scale ratio-plot to error of the reference histogram (1 standard deviation)

    Returns
    -------
    dict[str, dict]
        A dict containing all dicts that are usually written to the dat file. The key is the analysis ID.

    Raises
    ------
    IOError
        If the program does not have read access to .plot or .yoda files, or if it cannot write the output .dat files.

    Notes
    -----
    TODO The keys in the returned dict always includes a / rather than being the actual output file name.
    The keys will therefore differ from the actual output file names when hier_output == False.
    To get the actual file names, / should be replaced by _ when hier_output == False.
    """

    if verbose:
        logging.basicConfig(level=logging.DEBUG)

    # TODO: more elegant solution for getting rc_params by refactoring _parse_args.
    #  Then the 4 lines below can be replaced by 1 line
    stylename, _, rc_params_dict = _parse_args([style])
    stylename = stylename[0]    # Convert list to str
    rc_params_dict = rc_params_dict[stylename]  # Convert dict of dicts to dict
    del rc_params_dict['Title']

    ## Add pwd to search paths
    if path_pwd:
        rivet.addAnalysisLibPath(os.path.abspath("."))
        rivet.addAnalysisDataPath(os.path.abspath("."))
    for path in analysispaths:
        rivet.addAnalysisLibPath(os.path.abspath(path))
        rivet.addAnalysisDataPath(os.path.abspath(path))

    # Split the input file names and the associated plotting options given on the command line into two separate lists
    filelist, filenames, plotoptions = _parse_args(args)

    ## Check that the files exist
    for f in filelist:
        if not os.access(f, os.R_OK):
            raise IOError("Error: cannot read from %s" % f)

    plotdirs = plotinfodirs
    plotdirs += [os.path.abspath(os.path.dirname(f)) for f in filelist]
    plotdirs += (rivet.getAnalysisPlotPaths() if rivetplotpaths else [])

    # Create a list of all histograms to be plotted, and identify if they are 2D histos (which need special plotting)
    refhistos, mchistos = _get_histos(filelist, filenames, plotoptions, path_patterns, path_unpatterns)

    hpaths = []
    for aos in mchistos.values():
        for p in aos.keys():
            ps = rivet.stripOptions(p)
            if ps and ps not in hpaths:
                hpaths.append(ps)

    # Unique list of analyses
    anas = list(set([x.split("/")[1] for x in hpaths]))

    ## Take reference data from the Rivet search paths, if there is not already
    if rivetrefs:
        refhistos2 = _get_rivet_ref_data(anas, path_patterns, path_unpatterns)
        refhistos2.update(refhistos)
        refhistos = refhistos2
    ## Purge unmatched ref data entries to save memory
    keylist = list(refhistos.keys())
    for refhpath in keylist:
        if refhpath not in hpaths:
            del refhistos[refhpath]

    # Write each file
    plot_info_dicts = {}
    for plot_id in hpaths:
        outputdict = _make_output(
            plot_id, plotdirs, config_files,
            mchistos, refhistos, reftitle,
            plotoptions, stylename, rc_params_dict, mc_errs,
            nRatioTicks, skipWeights, removeOptions, deviation, 
            canvasText, verbose
        )
        if 'histograms' in outputdict: # protection against Counters
            plot_info_dicts[plot_id] = outputdict

    return anas, plot_info_dicts

