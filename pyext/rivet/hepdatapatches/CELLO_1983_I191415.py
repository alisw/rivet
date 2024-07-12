import yoda,math

def patch(path, ao):
    if ( "CELLO_1983_I191415" in path and
         ("d01" in path or "d02" in path or "d03" in path or
          "d04" in path or "d05" in path or "d06" in path)) :
        # remove the extra point in the overlap region
        if("d04" in path or "d05" in path or "d06" in path) :
            newAO = yoda.core.Scatter2D()
            newAO.setPath(ao.path())
            for i in range(0,len(ao.points())) :
                if( (i==5 and "d04" in path) or
                    ((i==5 or i==7) and "d04" not in path )) :
                    y = 0.5*(ao.points()[i].y()+ao.points()[i+1].y())
                    ylow = math.sqrt(ao.points()[i  ].yErrs()[0]**2+\
                                     ao.points()[i+1].yErrs()[0]**2)
                    yupp = math.sqrt(ao.points()[i  ].yErrs()[1]**2+\
                                     ao.points()[i+1].yErrs()[1]**2)
                    newAO.addPoint(ao.points()[i].x(),y,(0.,0.),(ylow,yupp))
                elif( (i==6 and "d04" in path) or
                      (i==6 or i==8) and "d04" not in path) :
                    continue
                else :
                    newAO.addPoint(ao.points()[i])
            ao=newAO
        # sort out the bin widths
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
            if("d01" in path or "d02" in path or "d03" in path) :
                if(i==6 or i==11 or i==15) :
                    xupp=xlow
                    temp=xupp
                elif(i==7 or i==12 or i==16) :
                    xlow = ao.points()[i].x()-ao.points()[i-1].x()-temp
                if("d01" in path) :
                    if(i==16) :
                        xupp=xlow
                else :
                    if(i==17) :
                        xupp=xlow
            elif("d04" in path) :
                if(i==0) :
                    xupp = 2.*xupp-4.35e-2
                    xlow=xupp
                    temp=xupp
                elif(i==1) :
                    xlow=4.35e-2
                elif(i==4) :
                    xupp = xlow
                    temp=xupp
                elif(i==5) :
                    xlow = 2.*xlow-temp
            elif("d05" in path) :
                if(i==0) :
                    xupp = ao.points()[i+1].x()-ao.points()[i].x()-2.75e-2
                    xlow = xupp
                elif(i==1) :
                    xlow=2.75e-2
                elif(i==3) :
                    xupp=xlow
                    temp=xupp
                elif(i==4) :
                    xlow = ao.points()[i].x()-ao.points()[i-1].x()-temp
                    xupp = xlow
                    temp= xupp
                elif(i==5) :
                    xlow=ao.points()[i].x()-ao.points()[i-1].x()-temp
                    xupp=xlow
                    temp=xlow
                elif(i==6) :
                    xlow= ao.points()[i].x()-ao.points()[i-1].x()-temp
            elif("d06" in path) :
                if(i==0) :
                    xupp = 2.*xupp-1.8e-2
                    xlow=xupp
                elif(i==1) :
                    xlow=1.8e-2
                elif(i==3) :
                    xupp=xlow
                    temp=xlow
                elif(i==4) :
                    xlow=2.*xlow-temp
                    xupp=xlow
                    temp=xupp
                elif(i==5) :
                    xlow=2.*xlow-temp
                    xupp=xlow
                    temp=xlow
                elif(i==6) :
                    xlow= 2*xlow-temp
                    xupp=2*xupp-1.3e-1
                elif(i==7) :
                    xlow=xupp
            ao.points()[i].setXErrs((xlow,xupp))


    return ao
