def patch(path, ao):
    if "MARKII_1988_I261194" in path and "d01" in path:
        step = 0.05
        for p in ao.points():
            if p.x()<0.45 :
                p.setXErrs(step)
                step=0.025
            elif p.x()>0.5 :
                p.setXErrs(step)
                step=0.05
            else :
                p.setXErrs((step,0.05))
                step=0.05
    return ao
