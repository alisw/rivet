import yoda
def patch(path, ao):
    # fix bin widths
    if "TOPAZ_1995_I381900" in path:
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
            if("d02" in path) :
                if(ao.points()[i].x()<2.23) :
                    xlow=0.15
                    xupp=0.15
                elif(ao.points()[i].x()>2.45 and ao.points()[i].x()<2.48) :
                    xlow=0.10
                    xupp=0.075
                elif(ao.points()[i].x()>2.60 and ao.points()[i].x()<3.22) :
                    xlow=0.075
                    xupp=0.075
                elif(ao.points()[i].x()>3.40 and ao.points()[i].x()<3.62) :
                    xlow=0.05
                    xupp=0.05
                elif(ao.points()[i].x()>3.76 and ao.points()[i].x()<3.78) :
                    xlow=0.11
                    xupp=0.11
                elif(ao.points()[i].x()>3.95 and ao.points()[i].x()<3.97) :
                    xlow=0.08
                    xupp=0.08
                elif(ao.points()[i].x()>4.13 and ao.points()[i].x()<4.15) :
                    xlow=0.10
                    xupp=0.10
                elif(ao.points()[i].x()>4.41 and ao.points()[i].x()<4.43) :
                    xlow=0.18
                    xupp=0.18
                elif(ao.points()[i].x()>4.70 and ao.points()[i].x()<4.72) :
                    xlow=0.11
                    xupp=0.08
                elif(ao.points()[i].x()>4.81 and ao.points()[i].x()<4.84) :
                    xlow=0.04
                    xupp=0.04
            ao.points()[i].setXErrs((xlow,xupp))
    return ao
