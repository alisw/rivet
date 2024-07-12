def patch(path, ao):
    # fix bin values, from bin number to cms energy
    x = [982.5,1016.,1056.,1098.]
    if "E605_1991_I302822" in path and ("d17" in path or "d18" in path or "d19" in path or
                                        "d20" in path or "d21" in path or "d22" in path):
        for p in ao.points() :
            p.setXErrs(0.1)
    return ao
