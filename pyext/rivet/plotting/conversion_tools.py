# Convert the old plot file format to a dict that can be parsed by plot2yaml.py
# Based on pyext/rivet/plotinfo.py
from __future__ import print_function
import os, re, logging
from rivet.util import texpand
import ast

pat_begin_block = re.compile(r'^(#*\s*)?BEGIN (\w+) ?(\S+)?')
pat_begin_name_block = re.compile(r'^(#*\s*)?BEGIN (\w+) ?(\S+)? ?(\w+)?')
pat_end_block =   re.compile(r'^(#*\s*)?END (\w+)')
pat_comment = re.compile(r'^\s*#|^\s*$')
pat_property = re.compile(r'^(\w+?)\s*=\s*(.*)$')
pat_property_opt = re.compile(r'^ReplaceOption\[(\w+=\w+)\]=(.*)$')
pat_path_property  = re.compile(r'^(\S+?)::(\w+?)=(.*)$')


def _is_end_marker(line, blockname):
    m = pat_end_block.match(line)
    return m and m.group(2) == blockname


def _is_comment(line):
    return pat_comment.match(line) is not None


def type_conversion(value):
    """Convert the value of a property into the correct type.

    Parameters
    ----------
    value : str
        Value of a property, retrieved from a .plot file

    Returns
    -------
    str | byte | number | tuple | list | dict | set | bool | None
        The input but converted to a different type.

    Raises
    ------
    TypeError
        If value is not of type str.
    """
    if not isinstance(value, str):
        raise TypeError('Expected value to be of type str but got type {}'.format(type(value)))

    try:
        # TODO: probably does not work for space-separated lists, such as "XCustomMajorTicks"
        return ast.literal_eval(value)
    except (ValueError, SyntaxError): # literal_eval might raise more errors
        return value


def convert_legacy_plotfile(filename, hpath, section='PLOT'):
    """Parse a plot file of the old format and return a dict of the same format as the new plot files.

    Parameters
    ----------
    filename : str
        Filename of the .plot file. Should include the .plot extension.
    hpath : str
        The histogram path, usually with the format /AnalysisID/HistogramID.
    section : str
        Either PLOT or HISTOGRAM. Mostly just PLOT
    Returns
    -------
    plot_settings : dict
        Dict with of the same format and style as the new .plot files.
        Returns an empty dict if nothing was found.

    Raises
    ------
    ValueError
        If section is not PLOT or HISTOGRAM.
        This includes the SPECIAL and FUNCTION sections, which could be parsed before.
    """
    if section not in ('PLOT', 'HISTOGRAM'):
        raise ValueError('Expected section to be PLOT or HISTOGRAM but got {}.'.format(section))

    ## Assemble the list of headers from any matching plotinfo paths and additional style files
    pat_paths = {}
    ret = {}

    if not os.access(filename, os.R_OK):
        return {}
    startreading = False
    with open(filename) as f:
        msec = None
        for line in f:
            m = pat_begin_block.match(line)
            if m:
                tag, pathpat = m.group(2,3)
                # pathpat could be a regex
                if pathpat not in pat_paths:
                    try:
                        pat_paths[pathpat] = re.compile(pathpat)
                    except TypeError:
                        logging.debug("Error reading plot file for {}. Skipping.".format(filename))
                        return {}
                if tag == section:
                    m2 = pat_paths[pathpat].match(hpath)
                    if m2:
                        msec = m2
                        startreading = True
                        continue
            if not startreading:
                continue
            if _is_end_marker(line, section):
                startreading = False
                continue
            elif _is_comment(line):
                continue
            vm = pat_property.match(line)

            if vm:
                prop, value = vm.group(1,2)
                if msec:
                    oldval = value
                    try:
                        ## First escape backslashes *except* regex groups, then expand regex groups from path match
                        value = value.encode("string-escape")
                        value = re.sub("(\\\\)(\\d)", "\\2", value) #< r-strings actually made this harder, since the \) is still treated as an escape!
                        value = msec.expand(value)
                    except Exception as e: # TODO: bad exception handling
                        value = oldval #< roll back escapes if it goes wrong
                ret[prop] = type_conversion(texpand(value)) #< expand TeX shorthands and convert type if necessary
            vm = pat_property_opt.match(line)
            if vm:
                prop, value = vm.group(1,2)
                ret['ReplaceOption[' + prop + ']'] = type_conversion(texpand(value))
    return ret
