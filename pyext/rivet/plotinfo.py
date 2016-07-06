import os, re
from .util import texpand

class PlotParser(object):
#class PlotStyler(object) or class PlotInfo(object):
    """Reads Rivet's .plot files and determines which attributes to apply to each histo path."""

    pat_begin_block = re.compile('^#*\s*BEGIN ([A-Z0-9_]+) ?(\S+)?')
    pat_end_block =   re.compile('^#*\s*END ([A-Z0-9_]+)')
    pat_comment = re.compile('^#|^\s*$')
    pat_property = re.compile('^(\w+?)=(.*)$')
    pat_path_property  = re.compile('^(\S+?)::(\w+?)=(.*)$')
    pat_paths = {}

    def __init__(self, plotpaths=None, addfiles=[]):
        """
        Parameters
        ----------
        plotpaths : list of str, optional
            The directories to search for .plot files.
            The default is to call the rivet.getAnalysisPlotPaths() function to get
            the directory where the .plot files can be found. (Usually equivalent to calling :command:`rivet-config --datadir`)

        Raises
        ------
        ValueError
            If `plotpaths` is not specified and calling
            :command:`rivet-config` fails.
        """
        self.addfiles = addfiles

        self.plotpaths = plotpaths
        if not self.plotpaths:
            try:
                import rivet
                self.plotpaths = rivet.getAnalysisPlotPaths()
            except Exception, e:
                sys.stderr.write("Failed to load Rivet analysis plot paths: %s\n" % e)
                raise ValueError("No plot paths given and the rivet module could not be loaded!")


    def getSection(self, section, hpath):
        """Get a section for a histogram from a .plot file.

        Parameters
        ----------
        section : ('PLOT'|'SPECIAL'|'HISTOGRAM')
            The section that should be extracted.
        hpath : str
            The histogram path, i.e. /AnalysisID/HistogramID .

        Todo
        ----
        Caching!
            At the moment the result of the lookup is not cached so every
            call requires a file to be searched for and opened.
        """
        if section not in ['PLOT', 'SPECIAL', 'HISTOGRAM']:
            raise ValueError("Can't parse section \'%s\'" % section)

        ## Decompose the histo path and remove the /REF prefix if necessary
        parts = hpath.strip('/').split('/')
        #if len(parts[0]) == 0:
        #    del parts[0]
        if parts[0] == "REF":
            del parts[0]
        if not parts:
            raise ValueError("Found empty histo path (or equal to /REF). Shouldn't be posible...")
        # if len(parts) == 1:
        #     parts.insert(0, "ANALYSIS")
        hpath = "/" + "/".join(parts[-2:])

        ## Assemble the list of headers from any matching plotinfo paths and additional style files
        base = parts[0] + ".plot"
        ret = {'PLOT': {}, 'SPECIAL': None, 'HISTOGRAM': {}}
        for pidir in self.plotpaths:
            plotfile = os.path.join(pidir, base)
            #print plotfile
            self._readHeadersFromFile(plotfile, ret, section, hpath)
            ## Don't break here: we can collect settings from multiple .plot files
            # TODO: So the *last* path wins? Hmm... reverse the loop order?
        # TODO: Also, is it good that the user-specific extra files override the official ones? Depends on the point of the extra files...
        for extrafile in self.addfiles:
            self._readHeadersFromFile(extrafile, ret, section, hpath)
        return ret[section]


    def _readHeadersFromFile(self, plotfile, ret, section, hpath):
        """Get a section for a histogram from a .plot file."""
        if not os.access(plotfile, os.R_OK):
            return
        startreading = False
        f = open(plotfile)
        for line in f:
            m = self.pat_begin_block.match(line)
            if m:
                tag, pathpat = m.group(1,2)
                #print tag, pathpat
                # pathpat could be a regex
                if not self.pat_paths.has_key(pathpat):
                    self.pat_paths[pathpat] = re.compile(pathpat)
                if tag == section:
                    if self.pat_paths[pathpat].match(hpath):
                        startreading = True
                        if section in ['SPECIAL']:
                            ret[section] = ''
                        continue
            if not startreading:
                continue
            if self.isEndMarker(line, section):
                startreading = False
                continue
            elif self.isComment(line):
                continue
            if section in ['PLOT', 'HISTOGRAM']:
                vm = self.pat_property.match(line)
                if vm:
                    prop, value = vm.group(1,2)
                    #print prop, value
                    ret[section][prop] = texpand(value)
            elif section in ['SPECIAL']:
                ret[section] += line
        f.close()


    def getHeaders(self, hpath):
        """Get the plot headers for histogram hpath.

        This returns the PLOT section.

        Parameters
        ----------
        hpath : str
            The histogram path, i.e. /AnalysisID/HistogramID .

        Returns
        -------
        plot_section : dict
            The dictionary usually contains the 'Title', 'XLabel' and
            'YLabel' properties of the respective plot.

        See also
        --------
        :meth:`getSection`
        """
        return self.getSection('PLOT', hpath)


    def getSpecial(self, hpath):
        """Get a SPECIAL section for histogram hpath.

        The SPECIAL section is only available in a few analyses.

        Parameters
        ----------
        hpath : str
            Histogram path. Must have the form /AnalysisID/HistogramID .

        See also
        --------
        :meth:`getSection`
        """
        return self.getSection('SPECIAL', hpath)


    def getHistogramOptions(self, hpath):
        """Get a HISTOGRAM section for histogram hpath.

        The HISTOGRAM section is only available in a few analyses.

        Parameters
        ----------
        hpath : str
            Histogram path. Must have the form /AnalysisID/HistogramID .

        See also
        --------
        :meth:`getSection`
        """
        return self.getSection('HISTOGRAM', hpath)


    def isEndMarker(self, line, blockname):
        m = self.pat_end_block.match(line)
        return m and m.group(1) == blockname


    def isComment(self, line):
        return self.pat_comment.match(line) is not None


    def updateHistoHeaders(self, hist):
        headers = self.getHeaders(hist.histopath)
        if headers.has_key("Title"):
            hist.title = headers["Title"]
        if headers.has_key("XLabel"):
            hist.xlabel = headers["XLabel"]
        if headers.has_key("YLabel"):
            hist.ylabel = headers["YLabel"]


def mkStdPlotParser(dirs=None, addfiles=[]):
    """
    Make a PlotParser with the standard Rivet .plot locations automatically added to
    the manually set plot info dirs and additional files.
    """
    if dirs is None:
        dirs = []
    from .core import getAnalysisPlotPaths
    dirs += getAnalysisPlotPaths()
    seen = set()
    dirs = [d for d in dirs if d not in seen and not seen.add(d)]
    return PlotParser(dirs, addfiles)
