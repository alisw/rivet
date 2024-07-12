def patch(path, ao):
    # fix bin widths
    if "DELPHI_2011_I890503" in path and "d02" in path:
        for p in ao.points():
            p.setXErrs(0.5)
    return ao
