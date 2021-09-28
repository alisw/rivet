def patch(path, ao):
    # fix bin widths
    if path == "/REF/JADE_1984_I202785/d03-x01-y01" :
        for p in ao.points() :
            p.setXErrs(0.15)
    return ao
