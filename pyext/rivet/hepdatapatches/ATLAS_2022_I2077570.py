
def patch(path, ao):
    needs_patching = [
        '/REF/ATLAS_2022_I2077570/d03-x01-y01',
        '/REF/ATLAS_2022_I2077570/d04-x01-y01',
        '/REF/ATLAS_2022_I2077570/d09-x01-y01', 
        '/REF/ATLAS_2022_I2077570/d10-x01-y01',
        '/REF/ATLAS_2022_I2077570/d13-x01-y01',
    ]
    if path in needs_patching:
      for p in ao.points():
          xLo, xHi = p.xErrs()
          p.setXErrs((xLo + 0.5, xHi + 0.5))
    return ao

