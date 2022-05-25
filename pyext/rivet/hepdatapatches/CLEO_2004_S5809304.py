def patch(path, ao):
    # fix bin widths
    if "CLEO_2004_S5809304" in path and "d01" in path :
        for p in ao.points():
            p.setXErrs(0.5)
    return ao
