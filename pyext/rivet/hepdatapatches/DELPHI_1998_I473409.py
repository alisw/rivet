def patch(path, ao):
    # fix bin widths
    if "DELPHI_1998_I473409" in path and ("d01" in path or "d02" in path or "d03" in path):
        for p in ao.points():
            p.setXErrs(0.5)
    return ao

