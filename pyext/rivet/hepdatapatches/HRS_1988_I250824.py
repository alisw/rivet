def patch(path, ao):
    # fix bin widths
    if path == "/REF/HRS_1988_I250824/d01-x01-y01" :
        for i in range(0,len(ao.points())) :
            if(i!=len(ao.points())-1) :
                xupp = 0.5*(ao.points()[i+1].x()-ao.points()[i].x())
            if(i!=0) :
                xlow = 0.5*(ao.points()[i].x()-ao.points()[i-1].x())
            if(i==0) :
                xlow=xupp
            elif(i==len(ao.points())-1) :
                xlow=xupp
            # hacks
            if(i==9) :
                xupp=xlow
                temp=xupp
            elif(i==10) :
                xlow = ao.points()[i].x()-ao.points()[i-1].x()-temp
            ao.points()[i].setXErrs((xlow,xupp))
    return ao
