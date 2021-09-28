def patch(path, ao):
    # fix bin widths and centroid
    if "BABAR_2007_S6895344" in path and ("d02" in path or "d04" in path):
        for p in ao.points():
            if "d02" in path : p.setX(10.54)
            p.setXErrs(0.5)
    return ao
