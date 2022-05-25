
def patch(path, ao):
    needs_patching = [ 
      '/REF/CLEO_1984_I193577/d01-x01-y01'
    ]
    if path in needs_patching:
      for p in ao.points():
          p.setErrs(1, (0.001, 0.001))
    return ao

