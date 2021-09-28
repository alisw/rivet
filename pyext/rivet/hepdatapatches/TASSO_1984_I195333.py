import yoda,math
# fix bin widths 
def patch(path, ao):
    if "TASSO_1984_I195333" in path and "d03" in path:
        for p in ao.points() :
            p.setXErrs(1)
    return ao
