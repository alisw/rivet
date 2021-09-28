
def patch(path, ao):
    needs_patching = [ 
      '/REF/UA5_1987_S1640666/d01-x01-y01',
      '/REF/UA5_1989_S1926373/d02-x01-y01',
      '/REF/UA5_1989_S1926373/d03-x01-y01',
    ]
    if path in needs_patching:
      bWidth = 0.5
      if 'd03' in path:
        bWidth = 1.0
      for p in ao.points():
          p.setErrs(1, (bWidth, bWidth))
    return ao

