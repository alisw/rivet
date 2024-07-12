def patch(path, ao):
    if "MARKII_1982_I177606" in path and ("d02" in path or "d03" in path):
        step = 0.1
        for p in ao.points():
            p.setXErrs(step)
    return ao
