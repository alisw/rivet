
def patch(path, ao):
    needs_patching = [ 
      '/REF/UA5_1982_S875503/d02-x01-y01',
      '/REF/UA5_1982_S875503/d02-x01-y02',
      '/REF/UA5_1982_S875503/d02-x01-y03',
      '/REF/UA5_1982_S875503/d02-x01-y04',
      '/REF/UA5_1982_S875503/d02-x01-y05',
      '/REF/UA5_1982_S875503/d02-x01-y06',
      '/REF/UA5_1982_S875503/d02-x01-y07',
      '/REF/UA5_1982_S875503/d02-x01-y08',
      '/REF/UA5_1982_S875503/d02-x01-y09',
      '/REF/UA5_1982_S875503/d02-x01-y10',
    ]
    if path in needs_patching:
      for p in ao.points():
          p.setErrs(1, (0.5, 0.5))
    return ao

