import yoda
def patch(path, ao):
    # remove average bins
    if "TASSO_1990_S2148048" in path :
        if ( "d02" in path or "d05" in path or
             "d06" in path or "d07" in path or "d08" in path ):
            newAO = yoda.core.Scatter2D()
            newAO.setPath(ao.path())
            for p in ao.points() :
                if("d02" in path and p.xMin()==0. and p.xMax()==1.) :
                    continue
                elif("d05" in path and p.xMin()==0. and p.xMax()==6.) :
                    continue
                elif("d06" in path and p.xMin()==0. and p.xMax()==1.) :
                    continue
                elif("d07" in path and p.xMin()==0. and p.xMax()==.3) :
                    continue
                elif("d08" in path and p.xMin()<0.7 and p.xMax()==1.) :
                    continue
                newAO.addPoint(p)
                ao=newAO
    return ao
