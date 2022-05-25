
def patch(path, ao):
    needs_patching = [
        '/REF/ATLAS_2020_I1790439/d23-x01-y01',
        '/REF/ATLAS_2020_I1790439/d25-x01-y01',
        '/REF/ATLAS_2020_I1790439/d28-x01-y01', 
        '/REF/ATLAS_2020_I1790439/d30-x01-y01',
        '/REF/ATLAS_2020_I1790439/d32-x01-y01',
        '/REF/ATLAS_2020_I1790439/d34-x01-y01',
        '/REF/ATLAS_2020_I1790439/d36-x01-y01',
        '/REF/ATLAS_2020_I1790439/d38-x01-y01',
        '/REF/ATLAS_2020_I1790439/d40-x01-y01',
        '/REF/ATLAS_2020_I1790439/d42-x01-y01',
        '/REF/ATLAS_2020_I1790439/d44-x01-y01',
    ]
    if path in needs_patching:
        if 'd23' in path:
            ao.point(0).setX(0.)
            ao.point(0).setXErrs(0.5)
            ao.point(1).setX(1.)
            ao.point(1).setXErrs(0.5)
            ao.point(2).setX(2.)
            ao.point(2).setXErrs(0.5)
            ao.point(3).setX(3.)
            ao.point(3).setXErrs(0.5)
        if 'd25' in path:
            ao.point(0).setX(0.)
            ao.point(0).setXErrs(0.5)
            ao.point(1).setX(1.)
            ao.point(1).setXErrs(0.5)
            ao.point(2).setX(2.)
            ao.point(2).setXErrs(0.5)
            ao.point(3).setX(3.)
            ao.point(3).setXErrs(0.5)
        if 'd28' in path or 'd30' in path:
            ao.point(0).setX(15.)
            ao.point(0).setXErrs(15.)
            ao.point(1).setX(45.)
            ao.point(1).setXErrs(15.)
            ao.point(2).setX(90.)
            ao.point(2).setXErrs(30.)
            ao.point(3).setX(235.)
            ao.point(3).setXErrs(115.)
        if 'd32' in path:
            ao.point(0).setX(-250.)
            ao.point(0).setXErrs(250.)
            ao.point(1).setX(60.)
            ao.point(1).setXErrs(60.)
            ao.point(2).setX(285.)
            ao.point(2).setXErrs(165.)
            ao.point(3).setX(1725.)
            ao.point(3).setXErrs(1275.)
        if 'd34' in path:
            ao.point(0).setX(-0.5)
            ao.point(0).setXErrs(0.5)
            ao.point(1).setX(0.5)
            ao.point(1).setXErrs(0.5)
            ao.point(2).setX(1.75)
            ao.point(2).setXErrs(0.75)
            ao.point(3).setX(5.75)
            ao.point(3).setXErrs(3.25)
        if 'd36' in path:
            ao.point(0).setX(-0.5)
            ao.point(0).setXErrs(0.5)
            ao.point(1).setX(0.7855)
            ao.point(1).setXErrs(0.7855)
            ao.point(2).setX(2.3565)
            ao.point(2).setXErrs(0.7855)
            ao.point(3).setX(3.927)
            ao.point(3).setXErrs(0.785)
            ao.point(4).setX(5.4975)
            ao.point(4).setXErrs(0.7855)
        if 'd38' in path or 'd40' in path:
            ao.point(0).setX(-30.)
            ao.point(0).setXErrs(30.)
            ao.point(1).setX(30.)
            ao.point(1).setXErrs(30.)
            ao.point(2).setX(90.)
            ao.point(2).setXErrs(30.)
            if 'd38' in path:
                ao.point(3).setX(235)
                ao.point(3).setXErrs(115)
        if 'd42' in path:
            ao.point(0).setX(60.)
            ao.point(0).setXErrs(60.)
            ao.point(1).setX(150.)
            ao.point(1).setXErrs(30.)
            ao.point(2).setX(200.)
            ao.point(2).setXErrs(20.)
            ao.point(3).setX(260.)
            ao.point(3).setXErrs(40.)
            ao.point(4).setX(350.)
            ao.point(4).setXErrs(50.)
            ao.point(5).setX(500.)
            ao.point(5).setXErrs(100.)
            ao.point(6).setX(1300.)
            ao.point(6).setXErrs(700.)
        if 'd44' in path:
            ao.point(0).setX(90.)
            ao.point(0).setXErrs(90.)
            ao.point(1).setX(250.)
            ao.point(1).setXErrs(70.)
            ao.point(2).setX(385.)
            ao.point(2).setXErrs(65.)
            ao.point(3).setX(525.)
            ao.point(3).setXErrs(75.)
            ao.point(4).setX(800.)
            ao.point(4).setXErrs(200.)
            ao.point(5).setX(1750.)
            ao.point(5).setXErrs(750.)
    return ao

