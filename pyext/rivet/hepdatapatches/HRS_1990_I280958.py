def patch(path, ao):
    # fix bin widths
    if "HRS_1990_I280958" in path and ("d03" in path or "d04" in path or
                                       "d05" in path or "d06" in path) :
        for i in range(0,len(ao.points())) :
            if(i!=len(ao.points())-1) :
                xupp = 0.5*(ao.points()[i+1].x()-ao.points()[i].x())
            if(i!=0) :
                xlow = 0.5*(ao.points()[i].x()-ao.points()[i-1].x())
            if(i==0) :
                xlow=xupp
            elif(i==len(ao.points())-1) :
                xupp=xlow
            # hacks
            if("d03" in path) :
                if(i==10) :
                    xupp=xlow
                    temp=xlow
                elif(i==11) :
                    xlow=2.*xlow-temp
            ao.points()[i].setXErrs((xlow,xupp))

    return ao
