def patch(path, ao):
    if "ALEPH_1996_I402895" in path and "d01" in path:
        for p in ao.points():
            p.setXErrs(0.5)
    return ao
