
def patch(path, ao):
    needs_patching = [ 
      '/REF/UA5_1989_S1926373/d01-x01-y01',
      '/REF/UA5_1989_S1926373/d02-x01-y01',
      '/REF/UA5_1989_S1926373/d03-x01-y01',
      '/REF/UA5_1989_S1926373/d04-x01-y01',
      '/REF/UA5_1989_S1926373/d05-x01-y01',
      '/REF/UA5_1989_S1926373/d06-x01-y01',
      '/REF/UA5_1989_S1926373/d07-x01-y01',
      '/REF/UA5_1989_S1926373/d08-x01-y01',
      '/REF/UA5_1989_S1926373/d09-x01-y01',
      '/REF/UA5_1989_S1926373/d10-x01-y01',
      '/REF/UA5_1989_S1926373/d11-x01-y01',
      '/REF/UA5_1989_S1926373/d11-x01-y02',
      '/REF/UA5_1989_S1926373/d12-x01-y01',
      '/REF/UA5_1989_S1926373/d12-x01-y02',
    ]
    if path in needs_patching:
      bWidth = 0.5
      if 'd01' in path or 'd02' in path:
        bWidth = 1.0
      for p in ao.points():
          p.setXErrs((bWidth, bWidth))
      if 'd01' in path:
        ao.point(28).setXErrs((1.0, 2.0))
        ao.point(30).setXErrs((6.0, 6.0))
      if 'd02' in path:
        ao.point(0).setXErrs((2.0, 2.0))
        ao.point(1).setXErrs((2.0, 1.0))
        ao.point(49).setXErrs((1.0, 2.0))
        ao.point(52).setXErrs((3.0, 3.0))
        ao.point(53).setXErrs((7.0, 7.0))
      if 'd03' in path:
        ao.point(12).setXErrs((0.5, 1.0))
        ao.point(13).setXErrs((3.0, 3.0))
      if 'd04' in path:
        ao.point(30).setXErrs((0.5, 1.0))
        ao.point(32).setXErrs((6.0, 6.0))
      if 'd05' in path:
        ao.point(48).setXErrs((0.5, 1.0))
        ao.point(51).setXErrs((1.0, 1.0))
        ao.point(52).setXErrs((9.0, 9.0))
      if 'd06' in path:
        ao.point(0).setXErrs((1.0, 1.5))
        ao.point(53).setXErrs((0.5, 1.0))
        ao.point(56).setXErrs((1.0, 1.0))
        ao.point(57).setXErrs((7.0, 7.0))
      if 'd07' in path:
        ao.point(20).setXErrs((0.5, 1.0))
        ao.point(21).setXErrs((1.0, 1.0))
        ao.point(22).setXErrs((3.5, 3.5))
      if 'd08' in path:
        ao.point(46).setXErrs((0.5, 1.0))
        ao.point(49).setXErrs((1.0, 1.0))
        ao.point(50).setXErrs((2.0, 2.0))
        ao.point(51).setXErrs((8.0, 8.0))
      if 'd09' in path:
        ao.point(75).setXErrs((0.5, 1.0))
        ao.point(79).setXErrs((1.0, 1.0))
        ao.point(80).setXErrs((1.5, 1.5))
        ao.point(81).setXErrs((3.0, 3.0))
        ao.point(82).setXErrs((11.5, 11.5))
      if 'd10' in path:
        ao.point(0).setXErrs((1.0, 1.0))
        ao.point(1).setXErrs((1.0, 0.5))
        ao.point(92).setXErrs((0.5, 1.0))
        ao.point(97).setXErrs((1.5, 1.5))
        ao.point(98).setXErrs((2.0, 2.0))
        ao.point(99).setXErrs((3.0, 3.0))
        ao.point(100).setXErrs((9.0, 9.0))
    return ao

