def patch(path, ao):
    # fix x values so not using bin number
    needs_patching = ["/REF/BESII_2004_I622224/d01-x01-y01",
                      "/REF/BESII_2004_I622224/d02-x01-y01",
                      "/REF/BESII_2004_I622224/d03-x01-y01",
                      "/REF/BESII_2004_I622224/d04-x01-y01",
                      "/REF/BESII_2004_I622224/d05-x01-y01",
                      "/REF/BESII_2004_I622224/d06-x01-y01",]
    if path in needs_patching:
        step=0.25
        for i in range(0,len(ao.points())):
            if(i==0) :
                ao.points()[i].setXErrs((0.25,0.2))
            else :
                ao.points()[i].setXErrs(step)
            if(i==0) : step = 0.075
            else     : step = 0.05
    return ao

