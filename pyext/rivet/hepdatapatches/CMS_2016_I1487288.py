
def patch(path, ao):
    needs_patching = [ 
      '/REF/CMS_2016_I1487288/d04-x01-y01'
    ]
    if path in needs_patching:
      for p in ao.points():
          xLo, xHi = p.xErrs()
          p.setErrs(1, (xLo + 0.5, xHi + 0.5))
    return ao

