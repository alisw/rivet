def patch(path, ao):
    # fix bin widths
    if "CLEOII_1999_I504672" in path and "d01" in path:
        for p in ao.points():
            p.setXErrs(0.5)
    return ao

