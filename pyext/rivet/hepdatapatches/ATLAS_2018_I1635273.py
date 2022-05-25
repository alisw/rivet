
def patch(path, ao):
    needs_patching = [ 
      '/REF/ATLAS_2018_I1635273/d01-x01-y01', 
      '/REF/ATLAS_2018_I1635273/d03-x01-y01'
      '/REF/ATLAS_2018_I1635273/d03-x01-y02'
      '/REF/ATLAS_2018_I1635273/d03-x01-y03'
      '/REF/ATLAS_2018_I1635273/d34-x01-y01'
    ]
    if path in needs_patching:
      for p in ao.points():
          p.setErrs(1, (0.5, 0.5))
    return ao

