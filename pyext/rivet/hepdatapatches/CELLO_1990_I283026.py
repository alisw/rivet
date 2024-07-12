def patch(path, ao):
    if ("CELLO_1990_I283026" in path and
        ("d01" in path or "d02" in path or "d03" in path)) :
        for i in range(0,len(ao.points())) :
            if(i!=len(ao.points())-1) :
                xupp = 0.5*(ao.points()[i+1].x()-ao.points()[i].x())
            if(i!=0) :
                xlow = 0.5*(ao.points()[i].x()-ao.points()[i-1].x())
            if(i==0) :
                xlow=xupp
            elif(i==len(ao.points())-1) :
                xupp=xlow
            if("d02" in path) :
                if(i==1) :
                    xupp=xlow
                    temp=xupp
                elif(i==2) :
                    xlow=2.*xlow-temp
                    xupp=xlow
                    temp=xupp
                elif(i==3) :
                    xlow=2.*xlow-temp
                    xupp=xlow
            elif("d03" in path) :
                if(i==1) :
                    xupp=xlow
                    temp=xupp
                elif(i==2) :
                    xlow=2.*xlow-temp
                    xupp=xlow
                    temp=xupp
                elif(i==3) :
                    xlow=2.*xlow-temp
                    xupp=xlow
                    temp=xupp
                elif(i==4) :
                    xlow=2.*xlow-temp
                    xupp=xlow
                    temp=xupp
                elif(i==5) :
                    xlow=2.*xlow-temp
                    xupp=xlow
                    temp=xupp
                elif(i==6) :
                    xlow=2.*xlow-temp
            elif("d01" in path) :
                if(i==0) :
                    xupp=2.*xupp-8.5e-3
                    xlow=xupp
                elif(i==1) :
                    xupp=2.*xupp-1.05e-2
                    xlow=xupp
                elif(i==2) :
                    xlow=xupp
                elif(i==4) :
                    xupp=xlow
                    temp=xlow
                elif(i==5) :
                    xlow=2.*xlow-temp
                    xupp=xlow
                    temp=xlow
                elif(i==6) :
                    xlow=2.*xlow-temp
                    xupp=xlow
                    temp=xlow
                elif(i==7) :
                    xlow=2.*xlow-temp
            ao.points()[i].setXErrs((xlow,xupp))
    return ao
