
def patch(path, ao):
    needs_patching = [ 
      '/REF/BELLE_2009_I823878/d01-x01-y01',
      '/REF/BELLE_2009_I823878/d01-x01-y02',
      '/REF/BELLE_2009_I823878/d01-x01-y03',
      '/REF/BELLE_2009_I823878/d01-x01-y04'
    ]
    if path in needs_patching:
      for p in ao.points():
          p.setErrs(1, (0.5, 0.5))
    return ao

