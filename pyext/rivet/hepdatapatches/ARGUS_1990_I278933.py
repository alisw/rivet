def patch(path, ao):
    # fix bin widths
    if "ARGUS_1990_I278933" in path and ("d01" in path or "d02" in path ):
        for p in ao.points():
            p.setXErrs(0.5)
    return ao
