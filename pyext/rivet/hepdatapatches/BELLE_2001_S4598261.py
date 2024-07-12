def patch(path, ao):
    # fix bin widths
    if "BELLE_2001_S4598261" in path and "d02" in path:
        for p in ao.points():
            p.setXErrs(0.5)
    return ao
