def patch(path, ao):
    # fix bin widths
    if "HRS_1985_I201482" in path and ("d19" in path or "d20" in path) :
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
            if("d20" in path) :
                if(i==1) :
                    xupp=2.*xupp-0.125
                    temp=xupp
                elif(i==2) :
                    xlow=0.125
                    xupp=0.125
                elif(i==3) :
                    xlow=2.*xlow-0.125
                elif(i==25) :
                    xlow=0.125
                    xupp=0.125
                elif(i==26) :
                    xlow=2.*xlow-0.125
            elif("d19" in path) :
                if(i==0) :
                    xupp=2.*xupp-0.25
                    xlow=xupp
                elif(i==1) :
                    xlow=0.25
                elif(i==2) :
                    xupp=2.*xupp-0.125
                elif(i==3) :
                    xlow=0.125
                elif(i==26) :
                    xupp=0.125
                elif(i==27) :
                    xlow=2.*xlow-.125
                elif(i==28) :
                    xupp=.25
                elif(i==29) :
                    xlow=2.*xlow-.25
                    xupp=xlow
            ao.points()[i].setXErrs((xlow,xupp))
    return ao
