import math
def patch(path, ao):
    # set bin widths
    if "AMY_1995_I406129" in path and "d01" not in path and "d05" not in path:
        for i in range(0,len(ao.points())) :
            if(i!=len(ao.points())-1) :
                xupp = 0.5*(ao.points()[i+1].x()-ao.points()[i].x())
            if(i!=0) :
                xlow = 0.5*(ao.points()[i].x()-ao.points()[i-1].x())
            if(i==0) :
                xlow=xupp
            elif(i==len(ao.points())-1) :
                xlow=xupp
            # # hacks
            if("d04" in path) :
                if(i==3 or i==10 or i==15) :
                    xupp=xlow
                    temp=xupp
                elif(i==4 or i==11 or i==16) :
                    xlow = 2.*xlow-temp
            elif("d02" in path or "d03" in path or "d06" in path) :
                if(i==7 or i==12) :
                    xupp=xlow
                    temp=xupp
                elif(i==8 or i==13) :
                    xlow = 2.*xlow-temp
            ao.points()[i].setXErrs((xlow,xupp))
    return ao
