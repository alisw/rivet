
def patch(path, ao):
    # fix bin widths
    if "CMS_2017_I1608166" in path:
        step=0.025
        for p in ao.points():
            p.setXErrs(step)
    return ao