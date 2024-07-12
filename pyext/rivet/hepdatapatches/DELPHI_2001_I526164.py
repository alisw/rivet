def patch(path, ao):
    if "DELPHI_2001_I526164" in path :
        # add bin width for plottings
        if "d04" in path :
            for p in ao.points() : p.setXErrs(0.5)
        # renumber bins as duplicate energy values
        elif "d03" in path :
            yval=30.
            if "y02" in path : yval=8.
            for p in ao.points() :
                p.setXErrs(0.5)
                if(p.x()<185) :
                    if(p.y()>yval) :
                        p.setX(1)
                    else :
                        p.setX(3)
                else :
                    if(p.y()>yval) :
                        p.setX(2)
                    else :
                        p.setX(4)
    return ao
