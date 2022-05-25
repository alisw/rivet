
def patch(path, ao):
    needs_patching = [ 
      '/REF/NMD_1974_I745/d01-x01-y01',
      '/REF/NMD_1974_I745/d01-x01-y02',
    ]
    if path in needs_patching:
      for p in ao.points():
          p.setErrs(1, (0.1, 0.1))
    return ao

