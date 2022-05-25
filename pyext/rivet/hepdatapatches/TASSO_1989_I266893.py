import yoda
# weird issue with bin low limits in hepdata, take upper limit of previous bin instead
def patch(path, ao):
    if path == "/REF/TASSO_1989_I266893/d12-x01-y01" :
        for i in range(1,len(ao.points())) :
            xmin = ao.points()[i-1].xMax()
            xmax = ao.points()[i].xMax()
            ao.points()[i].setX(0.5*(xmin+xmax))
            ao.points()[i].setXErrs(0.5*(xmax-xmin))
    return ao
