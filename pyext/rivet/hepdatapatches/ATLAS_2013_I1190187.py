def patch(path, ao):
    needs_patching = [u'/REF/ATLAS_2013_I1190187/d01-x01-y01',
                      u'/REF/ATLAS_2013_I1190187/d01-x01-y02',
                      u'/REF/ATLAS_2013_I1190187/d01-x01-y03',
                      u'/REF/ATLAS_2013_I1190187/d02-x01-y01',
                      u'/REF/ATLAS_2013_I1190187/d02-x01-y02',
                      u'/REF/ATLAS_2013_I1190187/d02-x01-y03',
                      u'/REF/ATLAS_2013_I1190187/d03-x01-y01']
    # set the x bin width (sum of the two x errors) to 1.0
    if path in needs_patching:
        for p in ao.points():
            p.setXErrs(0.5)
    return ao
