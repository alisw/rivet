def patch(path, ao):
    needs_patching = [u'/REF/LHCB_2016_I1454404/d04-x01-y01',
                      u'/REF/LHCB_2016_I1454404/d04-x01-y02',
                      u'/REF/LHCB_2016_I1454404/d05-x01-y01',
                      u'/REF/LHCB_2016_I1454404/d05-x01-y02',
                      u'/REF/LHCB_2016_I1454404/d06-x01-y01',
                      u'/REF/LHCB_2016_I1454404/d06-x01-y02',
                      u'/REF/LHCB_2016_I1454404/d07-x01-y01',
                      u'/REF/LHCB_2016_I1454404/d08-x01-y01',
                      u'/REF/LHCB_2016_I1454404/d09-x01-y01',
                      u'/REF/LHCB_2016_I1454404/d10-x01-y01']
    if path in needs_patching:
        for p in ao.points():
            bw = p.errPlus(1)+p.errMinus(1)
            p.scaleY(1./bw)
    return ao
