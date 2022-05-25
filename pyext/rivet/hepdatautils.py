from __future__ import print_function

def compare_with_hepdata(yodafile, inspire_id=0, yodafile_from_hepdata=None, output=None):
    """\
    Compare a YODA reference data file, intended for inclusion in Rivet, with the YODA file downloaded from HEPData.
    Make the comparison using the yodadiff script distributed with YODA (https://yoda.hepforge.org/trac/browser/bin/yodadiff).

    :param yodafile: name of YODA reference data file (intended for inclusion in Rivet)
    :param inspire_id: INSPIRE ID (to download the YODA file from HEPData)
    :param yodafile_from_hepdata: name of YODA file already downloaded from HEPData
    :param output: name of output file for yodadiff instead of stdout (note: -o option of yodadiff not implemented)
    :return: True or False depending on whether YODA files are compatible
    """

    if inspire_id:
        # Extract Rivet analysis name from last pathname component of yodafile (and discard extension).
        import os
        tail = os.path.split(yodafile)[1]
        rivet_analysis_name = os.path.splitext(tail)[0]
        yodafile_from_hepdata = download_from_hepdata(inspire_id, rivet_analysis_name)

    if yodafile_from_hepdata:
        print('Comparing {} with {}'.format(yodafile, yodafile_from_hepdata))
        import subprocess
        args = ['yodadiff', yodafile, yodafile_from_hepdata]
        if output:
            print('Writing output of "{}" to {}'.format(' '.join(args), output))
            with open(output, 'w') as out:
                # Note: -o|--output option of yodadiff is not implemented.
                exit_status = subprocess.call(args, stdout=out, stderr=out)
        else:
            exit_status = subprocess.call(args)
        if exit_status:
            return False
    else:
        print('No YODA file from HEPData specified!')
        return False

    return True


def download_from_hepdata(inspire_id, rivet_analysis_name=None, prefix='.'):
    """\
    Download the latest YODA reference data file from HEPData identified by the INSPIRE ID.
    Optionally pass the Rivet analysis name to write in the YODA file exported from HEPData.

    :param inspire_id: INSPIRE ID
    :param rivet_analysis_name: Rivet analysis name to override default
    :return: name of YODA file downloaded from HEPData
    """

    import tarfile, io, os, requests, gzip

    hdurl = 'https://hepdata.net/record/ins{}'.format(inspire_id)
    payload = {'format': 'yoda', 'rivet': rivet_analysis_name}
    response = requests.get(hdurl, params=payload)
    if response.history:
        print('Downloading from {}'.format(response.history[0].url))
    else:
        print('Downloading from {}'.format(response.url))

    if response.status_code != 200:
        print('Download failed ({} {}), does {} exist?'.format(response.status_code, response.reason, hdurl))
        return None

    try:
        with tarfile.open(mode='r:gz', fileobj=io.BytesIO(response.content)) as tar:
          #tar.extractall()
          #yodafile_from_hepdata = tar.getnames()[0]
          yodafile_from_hepdata = tar.members[0].name + '.gz'
          yodapath = prefix + '/' + yodafile_from_hepdata
          with gzip.open(yodapath, 'wb') as fout:
              fout.write(tar.extractfile(tar.members[0].name).read())
          os.chmod(yodapath, 0o644)
    except tarfile.TarError as e:
        print('Error reading tarfile ({})'.format(str(e)))
        return None

    print('Downloaded {}'.format(yodafile_from_hepdata))
    return yodapath



def patch_yodaref(yoda_from_hepdata, pattern=None, unpattern=None):
    """\
    Take a YODA file and check if the reference data contained in the file is in need of post-processing.
    If so, apply relevant post-processing steps and return.

    :param yoda_from_hepdata: YODA filename containing reference data from HEPData for post-processing
    :param pattern: optional positive-filtering regex to pass to yoda.read(). Empty str = None
    :param unpattern: optional negative-filtering regex to pass to yoda.read(). Empty str = None
    :return: dict of post-processed YODA data objects from HEPData
    """

    import importlib, yoda
    import rivet.hepdatapatches as hdpatch
    # Get analysis object, request Ref(Un)match strings to pass as pattern & unpattern
    hepdata_content = yoda.read(yoda_from_hepdata, True, pattern, unpattern)
    for tag in hepdata_content:
        if tag.startswith("/REF"):
            routine, tableid = tag.rstrip("/")[5:].split('/')
            if hasattr(hdpatch, routine):
                # get relevant patch function for this routine and apply patch
                routine_patcher = importlib.import_module("rivet.hepdatapatches." + routine)
                hepdata_content[tag] = routine_patcher.patch(tag, hepdata_content[tag])
            # The following hack is to remove the ErrorBreakdown for HepData entries
            # where the histogram has inconsistent breakdowns across bins
            # This is intended as a temporary workaround until we've fixed all the entries
            try:
              try:
                  nvars = len(hepdata_content[tag].variations())
              except RuntimeError:
                  raise AssertionError("WARNING - Table {0} for {1} has a broken ErrorBreakdown!".format(tableid, routine))
              for p in hepdata_content[tag].points():
                  assert(nvars == len(p.errMap().keys()))
            except AssertionError:
                print ("WARNING - Table {0} for {1} has an inconsistent ErrorBreakdown across points!".format(tableid, routine))
                hepdata_content[tag].rmAnnotation('ErrorBreakdown')
                hepdata_content[tag].rmVariations()
    return hepdata_content
