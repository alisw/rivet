
def patch(path, ao):
    needs_patching = [ 
      '/REF/TPC_1987_I235694/d02-x01-y03',
      '/REF/TPC_1987_I235694/d02-x01-y04',
      '/REF/TPC_1987_I235694/d03-x01-y01',
      '/REF/TPC_1987_I235694/d03-x01-y02',
      '/REF/TPC_1987_I235694/d03-x01-y03',
      '/REF/TPC_1987_I235694/d03-x01-y04',
      '/REF/TPC_1987_I235694/d04-x01-y01',
      '/REF/TPC_1987_I235694/d04-x01-y02',
      '/REF/TPC_1987_I235694/d04-x01-y03',
      '/REF/TPC_1987_I235694/d04-x01-y04',
      '/REF/TPC_1987_I235694/d05-x01-y01',
      '/REF/TPC_1987_I235694/d05-x01-y02',
      '/REF/TPC_1987_I235694/d05-x01-y03',
      '/REF/TPC_1987_I235694/d05-x01-y04',
    ]
    if path in needs_patching:
      for p in ao.points():
        p.setErrs(1, (0.5, 0.5))
    return ao

