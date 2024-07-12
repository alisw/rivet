def patch(path, ao):
    # fix bin widths
    if "HRS_1986_I18502" in path and ("d01" in path or "d02" in path or "d03" in path) :
        step=1.
        if("d03" in path) : step=0.5
        for p in ao.points() :
            p.setXErrs(step)
    return ao
