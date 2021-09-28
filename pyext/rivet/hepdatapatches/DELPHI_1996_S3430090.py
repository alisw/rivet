def patch(path, ao):
    # fix bin widths
    if "DELPHI_1996_S3430090" in path and ("d35" in path or "d36" in path or "d37" in path or
                                           "d38" in path or "d39" in path or "d40" in path):
        for p in ao.points():
            p.setXErrs(0.5)
    return ao

