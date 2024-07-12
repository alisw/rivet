from subprocess import Popen, PIPE
from sys import stdout, argv
from glob import glob

# TODO
# get rid of shell=True

def printINFO(msg):
    for line in msg.splitlines():
        print "INFO: %s" % line
        continue

    stdout.flush()

    return


def printERROR(msg):
    for line in msg.splitlines():
        print "ERROR: %s" % line
        continue

    stdout.flush()

    return


def allRivetAnalyses(env):
    cmd = ["rivet", "--list-analyses"]

    printINFO(" ".join(cmd))

    p = Popen(cmd,
            stdout=PIPE,
            env=env)

    return p.communicate()[0].split()


def getRivetEnv(rivetPath):
    # get environments
    cmd = [".", "%srivetenv.sh;" % rivetPath, "env", "-0"]
    p = Popen(" ".join(cmd),
            shell=True,
            stderr=PIPE,
            stdout=PIPE)

    pout = p.communicate()[0]
    env = {}
    for line in pout.split('\x00'):
        if not line: continue
        k, v = line.split('=', 1)
        env[k] = v

        continue

    return env


def checkAnalysis(oldRivetEnv, newRivetEnv, analysisName, hepmcgzFile):

    # first, run old analysis
    cmd = ["zcat", hepmcgzFile]
    printINFO(" ".join(cmd))

    pOldEvents = Popen(cmd,
            stdout=PIPE,
            stderr=PIPE)

    cmd = ["rivet", "-a", analysisName, "-H", "%s.old.yoda" % analysisName]
    printINFO(" ".join(cmd))

    logOld = open("%s.old.log" % analysisName, 'w')
    pOld = Popen(cmd,
            stdin=pOldEvents.stdout,
            stderr=logOld,
            stdout=logOld,
            env=oldRivetEnv)

    pOldEvents.stdout.close()


    # first, run new analysis
    cmd = ["zcat", hepmcgzFile]
    printINFO(" ".join(cmd))

    pNewEvents = Popen(cmd,
            stdout=PIPE,
            stderr=PIPE)

    cmd = ["rivet", "-a", analysisName, "-H", "%s.new.yoda" % analysisName]
    printINFO(" ".join(cmd))

    logNew = open("%s.new.log" % analysisName, 'w')
    pNew = Popen(cmd,
            stdin=pNewEvents.stdout,
            stderr=logNew,
            stdout=logNew,
            env=newRivetEnv)

    pNewEvents.stdout.close()


    # wait for both analyses to finish
    if pOld.wait() or pNew.wait():
        printERROR("analysis %s failed to run." % analysisName)
        return True

    cmd = ["yodadiff",
            "%s.old.yoda" % analysisName,
            "%s.new.yoda" % analysisName]

    printINFO(" ".join(cmd))

    # run yodadiff on resulting histograms
    pYodaDiff = Popen(cmd,
            stdout=PIPE,
            stderr=PIPE,
            env=newRivetEnv)

    output = pYodaDiff.communicate()[0]

    return bool(output)


def main():

    oldRivetEnv = getRivetEnv(argv[1].rstrip("/") + "/")
    newRivetEnv = getRivetEnv(argv[2].rstrip("/") + "/")

    allOldAnalyses = set(allRivetAnalyses(oldRivetEnv))
    allNewAnalyses = set(allRivetAnalyses(newRivetEnv))

    printINFO("--------------------------------")
    printINFO("--------------------------------")
    printINFO("old rivet executable and version:")
    p = Popen(["which", "rivet"], env=oldRivetEnv)
    p.wait()
    p = Popen(["rivet", "--version"], env=oldRivetEnv)
    p.wait()
    printINFO("--------------------------------")
    printINFO("new rivet executable and version:")
    p = Popen(["which", "rivet"], env=newRivetEnv)
    p.wait()
    p = Popen(["rivet", "--version"], env=newRivetEnv)
    p.wait()
    printINFO("--------------------------------")
    printINFO("--------------------------------")

    # read in combinations file and fill the analysis-to-testfile map
    hepmcMap = {}
    mapfile = open("/hepforge/home/rivet/hepmc-refdata/InFiles/COMBINATIONS", 'r')
    for line in mapfile:
        hepmcFile, analysis = line.split()[:2]

        hepmcMap[analysis] = "/hepforge/home/rivet/hepmc-refdata/%s" % hepmcFile
        continue

    for analysis, hepmcFile in hepmcMap.iteritems():
        printINFO("")
        printINFO("--------------------------------")
        printINFO("%s" % analysis)
        printINFO("--------------------------------")
        printINFO("datafile: %s" % hepmcFile)
        printINFO("")

        if analysis not in allOldAnalyses:
            printINFO("analysis %s not implemented in rivet at\n%s" % \
                    (analysis, argv[1]))
            continue

        if analysis not in allNewAnalyses:
            printINFO("analysis %s not implemented in rivet at\n%s" % \
                    (analysis, argv[2]))
            continue


        if checkAnalysis(oldRivetEnv, newRivetEnv, analysis, hepmcFile):
            printERROR("analysis %s fails histogram comparison." % analysis)
        else:
            printINFO("analysis %s passes histogram comparison." % analysis)

        continue

    return 0


if __name__ == "__main__":
    main()
