
def patch(path, ao):
    needs_patching = [ 
      '/REF/DM2_1990_I297706/d01-x01-y01',
      '/REF/DM2_1990_I297706/d02-x01-y01'
    ]
    if path in needs_patching:
      for p in ao.points():
          p.setErrs(1, (0.5, 0.5))
    return ao

