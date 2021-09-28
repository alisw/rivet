
def patch(path, ao):
    needs_patching = [ 
      '/REF/CRYSTAL_BALL_1988_I261078/d01-x01-y01'
    ]
    if path in needs_patching:
      for p in ao.points():
          p.setErrs(1, (0.1, 0.1))
    return ao

