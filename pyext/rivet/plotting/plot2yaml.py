from __future__ import print_function
import os, re, io, logging
import rivet, yoda
from rivet.plotting.conversion_tools import convert_legacy_plotfile

#def literal_presenter(dumper, data):
#    """Function that tells the yaml writer how to output yoda histograms into a yaml file.
#    Currently uses yoda.writeYODA
#
#    Parameters
#    ----------
#    dumper : yaml.YAML
#        The file writer
#    data : yoda histogram
#        The histogram that will be written to the file.
#
#    Returns
#    -------
#    str (?)
#        The str that will be written to the file
#    """
#    with io.StringIO() as filelike_str:
#        yoda.writeYODA(data, filelike_str)
#        output_str = filelike_str.getvalue()
#    return dumper.represent_scalar('tag:yaml.org,2002:str', output_str, style='|')
#
#yaml.add_multi_representer(yoda.AnalysisObject, literal_presenter)


def _parse_yaml_plotfile(filename, hpath):
    """Parse a .plot file with the yaml format and return the settings.

    Parameters
    ----------
    filename : str
        The name of a .plot file with yaml-like syntax.
    hpath : str
        The histogram path, usually with the format /AnalysisID/HistogramID .

    Returns
    -------
    plotfile_configs : dict
        All the settings from the file
    """
    plotfile_configs = {}
    with open(filename, 'r') as file:
        for configs in yoda.util._parseyaml(file):
            if re.match(configs['name'], hpath):   # TODO: old code uses match instead of fullmatch. Is this intentional?
                plotfile_configs.update(configs['plot features'])
    return plotfile_configs


def _is_legacy_format(filename):
    """Check if a file contains the old .plot format or the new yaml format.

    Parameters
    ----------
    filename : str
        Path to the .plot file that will be checked.

    Returns
    -------
    is_legacy : bool
        True if the file is the old format. False otherwise.

    Notes
    -----
    This does not make an explicit check that the format is correct yaml syntax.
    Instead, it only checks whether the file contains #BEGIN PLOT or not.
    """
    is_legacy = False
    with open(filename) as file:
        for line in file:
            if re.match(r'^(#*\s*)?BEGIN (\w+) ?(\S+)? ?(\w+)?', line):
                is_legacy = True
                break
    return is_legacy


def _get_matching_plot_configs_from_file(hpath, plotfilepath):
    """Open plotfilepath and get the settings with the corresponding hpath ID.

    Parameters
    ----------
    hpath : str
        The histogram path, usually with the format /AnalysisID/HistogramID.
    plotfilepath : str
        The path to a .plot file, either with yaml syntax or in the old format.

    Returns
    -------
    plot_configs : dict
        Dictionary of settings for the corresponding histogram from the .plot file.
        Will return an empty dict if the hpath could not be found in the file.
    """
    plot_configs = {}
    if not os.access(plotfilepath, os.R_OK):
        return {}

    if _is_legacy_format(plotfilepath):
        new_plot_settings = convert_legacy_plotfile(plotfilepath, hpath)
    else:
        new_plot_settings = _parse_yaml_plotfile(plotfilepath, hpath)

    if new_plot_settings:
        logging.debug('Found file {} with plot settings for histogram with ID {}'.format(plotfilepath, hpath))

    plot_configs.update(new_plot_settings)

    return plot_configs


def get_plot_configs(hpath, plotdirs=[], config_files=[]):
    """Get all settings for the hpath analysis by reading through the settings of plotdirs and config_files

    Parameters
    ----------
    hpath : str
        The histogram path,  usually with the format /AnalysisID/HistogramID .
    plotdirs : list[str]
        Directories containing .plot files. The settings in the files with name AnalysisID.plot will be parsed and added to plot_configs.
    config_files : Iterable[str]
        Extra .plot files with settings. If there are sections in these yaml files with name /AnalysisID/HistogramID, these settings will be added to plot_configs.

    Returns
    -------
    plot_configs : dict
        Dict containing all the settings from the .plot files that match the hpath ID.

    Note
    ----
    If a key exists in multiple configs, the last setting will included in plot_configs.
    As an example, is multiple .plot files have the `xlabel` setting, the xlabel specified in the last file will be used.
    """
    # Remove duplicates
    plotdirs = list(set(plotdirs))

    id_parts = rivet.AOPath(hpath).basepathparts()

    plot_configs = {}
    for plotdir in plotdirs:
        plotfilepath = os.path.join(plotdir, id_parts[0] + '.plot')
        plotfile_configs = _get_matching_plot_configs_from_file(hpath, plotfilepath)   # TODO: Can I pass hpath here or will that cause problems with /REF?
        plot_configs.update(plotfile_configs)

    for plotfilepath in config_files:
        plotfile_configs = _get_matching_plot_configs_from_file(hpath, plotfilepath)   # TODO: Can I pass hpath here or will that cause problems with /REF?
        plot_configs.update(plotfile_configs)

    return plot_configs

