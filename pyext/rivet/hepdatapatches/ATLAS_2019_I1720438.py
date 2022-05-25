
def patch(path, ao):
    needs_patching = [ 
      '/REF/ATLAS_2019_I1720438/d20-x01-y01',
    ]
    if path in needs_patching:
      for i in range(ao.numPoints()):
          ao.point(i).setVal(1, float(i))
          ao.point(i).setErrs(1, (0.5, 0.5))
    return ao

