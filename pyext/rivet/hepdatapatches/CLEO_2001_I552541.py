def patch(path, ao):
    # fix bin widths
    if "CLEO_2001_I552541" in path and ("d03" in path or "d04" in path):
        for p in ao.points():
            p.setXErrs(0.5)
    return ao
