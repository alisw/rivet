def patch(path, ao):
    # fix bin widths and centroid
    if "BABAR_2005_S6181155" in path and "d03" in path:
        for p in ao.points():
            p.setXErrs(0.5)
    return ao
