import yoda
def patch(path, ao):
    if( "TASSO_1988_I263859" in path) :
        newAO = yoda.core.Scatter2D()
        newAO.setPath(ao.path())
        for p in ao.points() :
            if("d01" in path and p.xMin()==0. and p.xMax()==1.) :
                continue
            elif("d02" in path and p.xMin()==0. and p.xMax()==.3) :
                continue
            elif("d03" in path and p.xMin()==0.5 and p.xMax()==1.) :
                continue
            elif("d04" in path and p.xMin()==0. and p.xMax()==2.5) :
                continue
            elif("d05" in path and p.xMin()==0. and p.xMax()==0.5) :
                continue
            elif("d06" in path and p.xMin()==2. and p.xMax()==34.) :
                continue
            elif("d07" in path and p.xMin()==0. and p.xMax()==5.) :
                continue
            elif("d08" in path and p.xMin()==0. and p.xMax()==2.) :
                continue
            elif("d09" in path and p.xMin()==0. and p.xMax()==5.) :
                continue
            elif("d10" in path and p.xMin()==0. and p.xMax()==1.) :
                continue
            elif("d11" in path and p.xMin()==0. and p.xMax()==5.) :
                continue
            if "d06" in path : p.setXErrs(1)
            newAO.addPoint(p)
        ao=newAO
    return ao
