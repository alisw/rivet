def patch(path, ao):
    # fix bin widths
    if "DELPHI_1999_I448370" in path and "d09" in path:
        for p in ao.points():
            p.setXErrs(0.5)
    return ao

