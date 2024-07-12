def patch(path, ao):
    # fix bin widths
    if "ARGUS_1991_I315059" in path and "d01" in path:
        for p in ao.points():
            p.setXErrs(0.5)
    return ao
