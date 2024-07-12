import yoda
def patch(path, ao):
    # fix x values, hepdata outputing bin number
    if "TPC_1988_I262143" in path and ("d06" in path or "d07" in path):
        if "d06" in path :
            bins = [0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,
                    0.1,0.11,0.12,0.13,0.14,0.16,0.18,0.2,0.22,0.25,0.3,0.35,0.4,0.5,0.6,0.7]
        else :
            bins = [0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,
                    0.1,0.11,0.25,0.3,0.35,0.4,0.5,0.6,0.7]
        for p in ao.points() :
            i = int(p.x())-1
            p.setX(0.5*(bins[i]+bins[i+1]))
            p.setXErrs(0.5*(bins[i+1]-bins[i]))
    return ao
