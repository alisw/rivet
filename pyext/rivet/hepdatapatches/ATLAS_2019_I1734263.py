def patch(path, ao):

    needs_patching = ['/REF/ATLAS_2019_I1734263/d01-x01-y01']

    if path in needs_patching:
        for p in ao.points():
            p.setErrs(1, (2.5, 2.5))
    return ao
