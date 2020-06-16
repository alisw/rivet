
def patch(path, ao):
    needs_patching = [ 
      '/REF/NA22_1986_I18431/d01-x01-y01',
      '/REF/NA22_1986_I18431/d02-x01-y01',
      '/REF/NA22_1986_I18431/d03-x01-y01',
    ]
    if path in needs_patching:
      import yoda
      newao = yoda.Scatter2D(ao.points()[1:])
      for anno in ao.annotations():
          newao.setAnnotation(anno, ao.annotation(anno))
      ao = newao
      for i in range(ao.numPoints()):
          ao.point(i).setErrs(1, (1.0, 1.0))
          if 'd01' in path:
              ao.point(i).setX(2*(i+1))
              ao.point(i).scaleY(1.0 / 20.94) # simga_inel(pi+ p) from DOI:10.1007/BF01550769
          elif 'd02' in path:
              ao.point(i).setX(2*(i+1))
              ao.point(i).scaleY(1.0 / 17.72) # simga_inel(K+ p) from DOI:10.1007/BF01550769
          elif 'd03' in path:
              ao.point(i).scaleY(1.0 / 32.40) # simga_inel(p p) from DOI:10.1007/BF01550769
    return ao

