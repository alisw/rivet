import yoda
def patch(path, ao):
    if( "CLEO_1985_I205668" in path and
        ("d01-x01" in path or "d02-x01" in path or "d03-x01" in path or "d04-x01" in path or "d05-x01" in path or
         "d06-x01" in path or "d07-x01" in path or "d08-x01" in path or "d09-x01" in path or "d10-x01" in path or
         "d11-x01" in path)) :
        # remove extra point
        if "d01-x01-y02" in path :
            newAO = yoda.core.Scatter2D()
            newAO.setPath(ao.path())
            for i in range(0,len(ao.points())) :
                if i==3 : continue
                newAO.addPoint(ao.points()[i])
            ao=newAO
        # bin widths
        for i in range(0,len(ao.points())) :
            if "d01-x01" in path :
                if(ao.points()[i].x()==.105) : continue
                if(ao.points()[i].x()<.2) :
                    ao.points()[i].setXErrs(0.01)
                elif(ao.points()[i].x()<.75) :
                    ao.points()[i].setXErrs(0.05)
                else :
                    ao.points()[i].setXErrs(0.1)
            elif  "d02-x01" in path:
                if(ao.points()[i].x()==.105) :
                    ao.points()[i].setXErrs(0.005)
                elif(ao.points()[i].x()!=0.06) :
                    ao.points()[i].setXErrs(0.01)
                else :
                    ao.points()[i].setXErrs(.06)
                    if "y01" in path :
                        ao.points()[i].setXErrs(0.03)
                    else :
                        ao.points()[i].setXErrs(0.04)
            elif  "d03-x01" in path:
                if(ao.points()[i].x()>0.16) :
                    ao.points()[i].setXErrs(0.015)
                elif(ao.points()[i].x()==.1475) :
                    ao.points()[i].setXErrs(.155-.1475)
                else:
                    ao.points()[i].setXErrs(.04)
            elif  "d04-x01" in path or "d07-x01" in path or "d08-x01" in path:
                ao.points()[i].setXErrs(.05)
            elif  "d05-x01" in path:
                if(ao.points()[i].x()!=.85) :
                    ao.points()[i].setXErrs(.025)
                else :
                    ao.points()[i].setXErrs(.05)
            elif  "d06-x01" in path:
                if(ao.points()[i].x()<.4) :
                    ao.points()[i].setXErrs(.025)
                elif(ao.points()[i].x()==.45) :
                    ao.points()[i].setXErrs(.05)
                elif(ao.points()[i].x()>=.575) :
                    ao.points()[i].setXErrs(.075)
            elif  "d09-x01" in path or "d10-x01" in path:
                if(ao.points()[i].x()>.17) :
                    ao.points()[i].setXErrs(.06)
                else :
                    ao.points()[i].setXErrs(.03)
            elif  "d11-x01" in path:
                if(ao.points()[i].x()==.29 or
                   ao.points()[i].x()==.48) :
                    ao.points()[i].setXErrs(.095)
                elif(ao.points()[i].x()==.76) :
                    ao.points()[i].setXErrs(.76-.575)
                elif(ao.points()[i].x()==.32) :
                    ao.points()[i].setXErrs(.04)
                elif(ao.points()[i].x()==.53) :
                    ao.points()[i].setXErrs(.17)
                else :
                    ao.points()[i].setXErrs(1.-.85)
    return ao
