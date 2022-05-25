def patch(path, ao):
    if "HRS_1986_I17781" in path and "d03" in path :
        # add bin width for plottings
            for p in ao.points() : p.setXErrs(0.5)
    return ao
