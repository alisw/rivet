def patch(path, ao):
    # fix bin widths
    if "L3_1992_I334954" in path :
        # set bin value and width for multiplicty dist and remove average bin
        if "d16" in path :
            x=6.
            for p in ao.points() :
                p.setXErrs(1.)
                p.setX(x)
                x+=2
            ao.rmPoint(len(ao.points())-1)
        # remove average bin from all
        # (its the one with the largest width as spans whole dist)
        else :
            iloc=-1
            eMax=0.
            for ix in range(0,len(ao.points())) :
                if(ao.points()[ix].xErrAvg()>eMax) :
                    eMax = ao.points()[ix].xErrAvg()
                    iloc=ix
            ao.rmPoint(iloc)
    return ao
