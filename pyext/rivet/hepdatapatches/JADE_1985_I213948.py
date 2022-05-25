import yoda
def patch(path, ao):
    # no bin widths, need to add them
    if "JADE_1985_I213948" in path :
        # photon spectra
        if "d01" in path or "d02" in path or "d03" in path:
            ii = [100,100,100,100,100,100]
            step = 0.1
            if "d01" in path :
                ii = [18,28,32,32,100,100]
                step=0.05
            elif "d02" in path:
                ii = [8,13,14,15,16,17]
            elif "d03" in path :
                ii = [7,12,13,14,15,16]
            for i in range(0,len(ao.points())) :
                if i<=ii[0]  :
                    ao.points()[i].setXErrs(0.005)
                elif i<=ii[1] :
                    ao.points()[i].setXErrs(0.01)
                elif i<=ii[2] :
                    ao.points()[i].setXErrs(0.025)
                elif i<=ii[3] :
                    ao.points()[i].setXErrs(0.075)
                elif i<=ii[4] :
                    ao.points()[i].setXErrs(step)
                elif i<=ii[5] :
                    ao.points()[i].setXErrs(0.2)
        elif "d04" in path or "d05" in path or "d06" in path or "d07" in path:
            if "d04" in path: 
                bins = [0.011,0.023,0.035,0.047,0.059,0.071,0.083,0.095,0.107,0.119,0.143,0.167,0.191,0.239]
            elif "d05" in path:
                bins = [0.018,0.036,0.054,0.072,0.108,0.144,0.216,0.360]
            elif "d06" in path:
                bins = [0.03,0.058,0.086,0.114,0.17,0.226,0.338,0.562]
            elif "d07" in path:
                bins = [0.04,0.056,0.088,0.28]
            for i in range(0,len(ao.points())) :
                ao.points()[i].setXErrs((ao.points()[i].x()-bins[i],bins[i+1]-ao.points()[i].x()))
    return ao
