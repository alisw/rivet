
def patch(path, ao):
    needs_patching = [ 
      '/REF/BABAR_2001_I558091/d01-x01-y01',
    ]
    if path in needs_patching:
      for p in ao.points():
          p.setErrs(1, (0.5, 0.5))
    return ao

