def patch(path, ao):
    if ( "ALEPH_1996_S3486095" in path and
         ("d44" in path or "d24" in path or
          "d23" in path or "d22" in path or
          "d21" in path or "d20" in path or
          "d19" in path or "d18" in path) ):
        step = 0.5
        if "d18" in path : step = 1.
        for p in ao.points():
            p.setXErrs(step)
    return ao
