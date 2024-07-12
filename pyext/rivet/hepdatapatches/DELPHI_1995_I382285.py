def patch(path, ao):
    # fix bin widths
    if "DELPHI_1995_I382285" in path and "d01" in path:
        for p in ao.points():
            p.setXErrs(0.5)
    return ao

