
def patch(path, ao):
    needs_patching = [
        '/REF/ATLAS_2018_I1711223/d12-x01-y01',
        '/REF/ATLAS_2018_I1711223/d20-x01-y01',
    ]
    if path in needs_patching:
        for i in range(ao.numPoints()):
            offset = 2.0 if 'd12' in path else 0.0
            ao.point(i).setX(float(i) + offset)
            ao.point(i).setXErrs(0.5)
    return ao

