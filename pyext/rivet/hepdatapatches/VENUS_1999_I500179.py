
def patch(path, ao):
    needs_patching = [ 
      '/REF/VENUS_1999_I500179/d01-x01-y01'
    ]
    if path in needs_patching:
      for p in ao.points():
          p.setErrs(1, (0.1, 0.1))
    return ao

