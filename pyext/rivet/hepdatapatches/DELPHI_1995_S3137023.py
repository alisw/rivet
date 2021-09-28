def patch(path, ao):
    # fix bin widths
    if "DELPHI_1995_S3137023" in path and ("d01" in path or "d04" in path or "d05" in path or
                                           "d06" in path or "d07" in path ):
        for p in ao.points():
            p.setXErrs(0.5)
    return ao

