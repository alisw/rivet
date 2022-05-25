def patch(path, ao):
    if "MARKII_1987_I247900" in path and "d01" not in path:
        step = 0.5
        for p in ao.points():
            p.setXErrs(step)
    return ao
