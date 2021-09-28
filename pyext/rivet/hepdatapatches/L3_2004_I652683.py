# bin widths
def patch(path, ao):
    if "L3_2004_I652683" in path:
        temp=path.split("/")[-1].split("-")
        i = int(temp[0][1:])
        if i==59:
            for p in ao.points():
                p.setXErrs(1)
        elif(i<=20) :
            for i in range(0,len(ao.points())) :
                if(i==0) :
                    ao.points()[i].setXErrs(5e-4)
                elif(i==1) :
                    ao.points()[i].setXErrs((5e-4,1e-3))
                elif(i<5) :
                    ao.points()[i].setXErrs(1e-3)
                elif(i==5) :
                    ao.points()[i].setXErrs((1e-3,5e-3))
                elif(i==6) :
                    ao.points()[i].setXErrs((5e-3,1e-2))
                else :
                    ao.points()[i].setXErrs(1e-2)
    return ao
