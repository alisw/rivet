
def patch(path, ao):
    needs_patching = [ 
      '/REF/ATLAS_2019_I1746286/d05-x01-y01',
      '/REF/ATLAS_2019_I1746286/d10-x01-y01',
      '/REF/ATLAS_2019_I1746286/d14-x01-y01',
      '/REF/ATLAS_2019_I1746286/d18-x01-y01',
    ]
    if path in needs_patching:
      for p in ao.points():
          xLo, xHi = p.xErrs()
          p.setErrs(1, (xLo + 0.5, xHi + 0.5))
    return ao

