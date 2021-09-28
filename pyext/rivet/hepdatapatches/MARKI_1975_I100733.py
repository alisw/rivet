# bin widths
def patch(path, ao):
    if "MARKI_1975_I100733" in path and "d03" in path :
        for p in ao.points() :
            if(p.x()==0.) :
                if "y01" in path :
                    if p.y()==0.9 :
                        p.setX(0.07)
                    else :
                        p.setX(0.09)
                elif "y02" in path :
                    if p.y()==3.8 :
                        p.setX(0.03)
                    elif p.y()==4.2 :
                        p.setX(0.05)
                    elif p.y()==5.54 :
                        p.setX(0.07)
                    else :
                        p.setX(0.09)
                elif "y03" in path :
                    if p.y()==13.1 :
                        p.setX(0.03)
                    elif p.y()==14.1 :
                        p.setX(0.05)
                    elif p.y()==14.8 :
                        p.setX(0.07)
                    else :
                        p.setX(0.09)
                p.setXErrs(0.01)
            elif (p.x()<0.40) :
                p.setXErrs(0.01)
            elif(p.x()<0.58) :
                p.setX(p.x()+0.01)
                p.setXErrs(0.02)
            else :
                if(p.x()!=0.84) : p.setX(p.x()+0.01)
                p.setXErrs(0.04)
    return ao
