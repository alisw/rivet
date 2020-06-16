
def patch(path, ao):
    needs_patching = [ 
      '/REF/BABAR_2006_I731865/d01-x01-y01',
      '/REF/BABAR_2006_I731865/d01-x01-y02'
    ]
    if path in needs_patching:
      for p in ao.points():
          p.setErrs(1, (0.5, 0.5))
    return ao

