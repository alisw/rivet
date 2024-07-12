import math
def patch(path, ao):
    # change so differenital in x not log10(x)
    if "ZEUS_2008_I780108" in path and "d20" in path:
        for p in ao.points():
            dxOld = math.log10(p.xMax())-math.log10(p.xMin())
            dxNew = p.xErrs()[0]+p.xErrs()[1]
            fact = dxOld/dxNew
            p.setY(p.y()*fact)
            for (key,value) in p.errMap().items() :
                temp=(value[0]*fact,value[1]*fact)
                p.setYErrs(temp,key)
    return ao
