def patch(path, ao):

    needs_patching = ['/REF/LHCB_2018_I1662483/d01-x01-y01']

    if path in needs_patching:
        for p in ao.points():
            p.setErrs(1, (0.5, 0.5))
    return ao
