import yoda
def patch(path, ao):
    # fix bin widths
    if ("OPAL_2004_S6132243" in path) :
        if("d15" in path or "d16" in path or "d17" in path or
           "d18" in path or "d19" in path or "d20" in path or
           "d21" in path or "d22" in path or "d23" in path or
           "d24" in path or "d25" in path or "d26" in path ) :
            for p in ao.points() :
                p.setXErrs(0.5)
        elif ("d27" in path) :
            ao.points()[0].setXErrs(39.)
    return ao
