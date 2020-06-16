
def patch(path, ao):
    needs_patching = [ 
      '/REF/MAC_1985_I206052/d01-x01-y01',
      '/REF/MAC_1985_I206052/d02-x01-y01',
      '/REF/MAC_1985_I206052/d02-x01-y02',
    ]
    if path in needs_patching:
      for p in ao.points():
          p.setErrs(1, (0.1, 0.1))
    return ao

