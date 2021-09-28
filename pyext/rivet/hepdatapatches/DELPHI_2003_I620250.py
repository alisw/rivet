def patch(path, ao):
    # fix bin widths
    if "DELPHI_2003_I620250" in path :
        if "d30" in path or "d31" in path:
            step=0.25
            for p in ao.points():
                p.setXErrs(step)
            ao.points()[-1].setXErrs(0.5)
        elif "d32" in path or "d33" in path:
            step=0.2
            for p in ao.points() : p.setXErrs(step)
        elif "d34" in path or "d35" in path:
            for p in ao.points() : 
                if(p.x()>.5 and p.x()<1.) :
                    p.setXErrs(0.125)
                elif(p.x()>.9 and p.x()<1.1) :
                    p.setXErrs(0.1)
                elif(p.x()>1.1 and p.x()<1.3) :
                    p.setXErrs(0.15)
                elif(p.x()>.2 and p.x()<.3) :
                    p.setXErrs(0.15)
                elif(p.x()<.2) :
                    p.setXErrs((p.x(),.1-p.x()))
                elif(p.x()>2 and p.x()<4) :
                    p.setXErrs(0.5)
                elif(p.x()>4 and p.x()<8) :
                    p.setXErrs(1)
                elif(p.x()>1.6 and p.x()<1.7) :
                    p.setXErrs((p.x()-1.4,2.-p.x()))
                else :
                    p.setXErrs(2)
        elif "d36" in path or "d37" in path:
            for p in ao.points() : 
                if(p.x()<0.6) :
                    p.setXErrs(0.1)
                elif(p.x()<0.8) :
                    p.setXErrs(0.125)
                elif(p.x()<1.1) :
                    p.setXErrs(0.175)
                elif(p.x()<1.9) :
                    p.setXErrs(0.2)
                else :
                    p.setXErrs(0.5)
        elif "d38" in path or "d39" in path:
            for p in ao.points() :
                if(p.x()<0.1) :
                    p.setXErrs(0.005)
                elif (p.x()<0.2) :
                    p.setXErrs(0.01)
                else :
                    p.setXErrs(0.02)
        elif ("d40" in path or "d41" in path or
              "d42" in path or "d43" in path):
            for p in ao.points() :
                if(p.x()<0.04) :
                    p.setXErrs(0.01)
                elif (p.x()<0.08) :
                    p.setXErrs(0.005)
                elif (p.x()<0.17) :
                    p.setXErrs(0.01)
                else :
                    p.setXErrs(0.02)
        elif "d44" in path or "d45" in path:
            for p in ao.points() :
                if(p.x()<0.2) :
                    p.setXErrs(0.01)
                else :
                    p.setXErrs(0.02)
        elif "d46" in path or "d47" in path:
            for p in ao.points() :
                if(p.x()<0.08) :
                    p.setXErrs(0.005)
                elif (p.x()<0.15) :
                    p.setXErrs(0.01)
                elif (p.x()<0.2) :
                    p.setXErrs(0.015)
                else :
                    p.setXErrs(0.02)
        elif "d48" in path or "d49" in path:
            for p in ao.points() :
                if(p.x()<0.11) :
                    p.setXErrs(0.005)
                elif(p.x()<0.21) :
                    p.setXErrs(0.01)
                else :
                    p.setXErrs(0.015)
        elif "d50" in path or "d51" in path:
            for p in ao.points() :
                if(p.x()<0.10) :
                    p.setXErrs(0.005)
                elif(p.x()<0.21) :
                    p.setXErrs(0.01)
                else :
                    p.setXErrs(0.02)
        elif "d52" in path or "d53" in path:
            for p in ao.points() :
                p.setXErrs(0.02)
        elif "d54" in path or "d55" in path:
            for p in ao.points() :
                if(p.x()<0.16) :
                    p.setXErrs(0.01)
                elif(p.x()<0.5) :
                    p.setXErrs(0.02)
                else :
                    p.setXErrs(0.03)
        elif ("d56" in path or "d57" in path or "d58" in path or
              "d59" in path or "d60" in path or "d61" in path or
              "d62" in path or "d63" in path ):
            for p in ao.points() :
                if(p.x()<0.06) :
                    p.setXErrs(0.005)
                elif(p.x()<0.16) :
                    p.setXErrs(0.01)
                else :
                    p.setXErrs(0.02)
        elif "d64" in path or "d65" in path:
            for p in ao.points() :
                if(p.x()<0.04) :
                    p.setXErrs(0.005)
                elif(p.x()<0.08) :
                    p.setXErrs(0.01)
                elif(p.x()<0.2) :
                    p.setXErrs(0.02)
                else :
                    p.setXErrs(0.025)
        elif ("d66" in path or "d67" in path or
              "d68" in path or "d69" in path ):
            for p in ao.points() :
                if(p.x()<0.06) :
                    p.setXErrs(0.005)
                elif(p.x()<0.12) :
                    p.setXErrs(0.01)
                elif(p.x()<0.2) :
                    p.setXErrs(0.02)
                elif(p.x()<0.41) :
                    p.setXErrs(0.025)
                else :
                    p.setXErrs(0.05)
        elif "d70" in path or "d71" in path:
            for p in ao.points() :
                if(p.x()>.04 and p.x()<.12) :
                    p.setXErrs(0.01)
                elif(p.x()>.12) :
                    p.setXErrs(0.02)
                elif(p.x()>0.02) :
                    p.setXErrs(0.005)
                elif(p.x()>0.017) :
                    p.setXErrs(0.002)
                elif(p.x()>0.006) :
                    p.setXErrs(0.003)
                else :
                    p.setXErrs(0.002)



    return ao
