import yoda,math
# removal of average bins
def patch(path, ao):
    # dists which need patching and the bin to remove
    needs_patching = { '/REF/AMY_1990_I283337/d02-x01-y01' : 5,
                       '/REF/AMY_1990_I283337/d03-x01-y01' : 5,
                       '/REF/AMY_1990_I283337/d04-x01-y01' : 5,
                       '/REF/AMY_1990_I283337/d05-x01-y01' : 4,
                       '/REF/AMY_1990_I283337/d06-x01-y01' : 5,
                       '/REF/AMY_1990_I283337/d07-x01-y01' : 3,
                       '/REF/AMY_1990_I283337/d08-x01-y01' : 5,
                       '/REF/AMY_1990_I283337/d09-x01-y01' : 4,
                       '/REF/AMY_1990_I283337/d10-x01-y01' : 8,
                       '/REF/AMY_1990_I283337/d11-x01-y01' : 8,
                       '/REF/AMY_1990_I283337/d12-x01-y01' : 4,
                       '/REF/AMY_1990_I283337/d14-x01-y01' : 4,
                       '/REF/AMY_1990_I283337/d15-x01-y01' : 4,
                       '/REF/AMY_1990_I283337/d16-x01-y01' : 4,
                       '/REF/AMY_1990_I283337/d17-x01-y01' : 3,
                       '/REF/AMY_1990_I283337/d18-x01-y01' : 3,
                       '/REF/AMY_1990_I283337/d19-x01-y01' : 4,
                       '/REF/AMY_1990_I283337/d20-x01-y01' : 2,
                       '/REF/AMY_1990_I283337/d21-x01-y01' : 5,
                       '/REF/AMY_1990_I283337/d22-x01-y01' : 5, }
    if path in needs_patching :
        newAO = yoda.core.Scatter2D()
        newAO.setPath(ao.path())
        ibin = needs_patching[path]
        for i in range(0,len(ao.points())) :
            if i !=ibin : newAO.addPoint(ao.points()[i])
        ao=newAO

    return ao
