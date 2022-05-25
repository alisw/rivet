def patch(path, ao):
    # fix bin widths
    if("HRS_1989_I276948") in path :
        if("d01" in path) :
            for i in range(0,len(ao.points())) :
                x = ao.points()[i].x()
                if(i>16) :
                    xmin = 0.1*(i-12)
                elif(i>0) :
                    xmin = 0.5*(ao.points()[i].x()+ao.points()[i-1].x())
                if(i>15) :
                    xmax = 0.1*(i-11)
                elif(i<len(ao.points())-1) :
                    xmax = 0.5*(ao.points()[i].x()+ao.points()[i+1].x())
                else :
                    xmax = 2.*ao.points()[i].x()-xmin
                if(i==0) : xmin=0.0725
                ao.points()[i].setXErrs((x-xmin,xmax-x))
        elif("d02" in path) :
            for i in range(0,len(ao.points())) :
                x = ao.points()[i].x()
                if(i==11) : xmin=0.7
                elif(i>8) : xmin=0.1*(i-4)
                elif(i>0) :
                    xmin = 0.5*(ao.points()[i].x()+ao.points()[i-1].x())
                if(i>9) : xmax=0.2*(i-7)+0.1
                elif(i>7) : xmax=0.1*(i-3)
                elif(i<len(ao.points())-1) :
                    xmax = 0.5*(ao.points()[i].x()+ao.points()[i+1].x())
                else :
                    xmax = 2.*ao.points()[i].x()-xmin
                if(i==0) : xmin=2.*ao.points()[i].x()-xmax
                ao.points()[i].setXErrs((x-xmin,xmax-x))
        else :
            for p in ao.points() :
                p.setXErrs(0.5)
    return ao
