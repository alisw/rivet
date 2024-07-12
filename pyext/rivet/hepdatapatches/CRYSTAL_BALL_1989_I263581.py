def patch(path, ao):
    # fix bin widths
    if "CRYSTAL_BALL_1989_I263581":
        for p in ao.points():
            p.setXErrs(0.025)
    return ao
