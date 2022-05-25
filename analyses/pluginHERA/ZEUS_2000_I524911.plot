BEGIN PLOT /ZEUS_2000_I524911/d01-x01-y01
LogY=0
Title=Differential $\phi$ distribution for $p_T = 0.5$ GeV
XLabel=$\phi (degrees)$ 
YLabel=$1/N dn/d\phi (rad^-1)$ 
XMin=-pi
XMax=pi
YMin=0.00
YMax=0.10
LegendYPos=0.4
LegendXPos=0.5

#YMajorTickMarks=5
#YMinorTickMarks=3

# + any additional plot settings you might like, see make-plots documentation
END PLOT

# ... add more histograms as you need them ...

BEGIN PLOT /ZEUS_2000_I524911/d01-x01-y02
LogY=0
Title=Differential $\phi$ distribution for $p_T = 1$ GeV
XLabel=$\phi (degree)$
YLabel=$1/N dn/d\phi (rad^-1)$
XMin=-pi
XMax=pi
YMin=0.00
YMax=0.08
YMajorTickMarks=10
YMinorTickMarks=3

# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /ZEUS_2000_I524911/d01-x01-y03
LogY=0
Title=Differential $\phi$ distribution for $p_T = 1.5$ GeV
XLabel=$\phi (degrees)$ 
YLabel=$1/N dn/d\phi (rad^-1)$ 
XMin=-pi
XMax=pi
YMin=0.00
YMax=0.06
YMajorTickMarks=5
YMinorTickMarks=3
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /ZEUS_2000_I524911/d01-x01-y04
Title=Differential $\phi$ distribution for $p_T = 2$ GeV
LogY=0
XLabel=$\phi (degrees)$ 
YLabel=$1/N dn/d\phi (rad^-1)$ 
XMin=-pi
XMax=pi
YMin=0.000
YMax=0.040
YMajorTickMarks=10
YMinorTickMarks=3
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /ZEUS_2000_I524911/d02-x01-y01
LogY=0
Title= Average $\cos \phi$ vrs $P_i$: 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0 GeV
XLabel=$P_i$ index
YLabel=$<\cos \phi>$
XMin=0
XMax=9.00
YMin=-0.15
YMax=0.10
YMajorTickMarks=10

# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /ZEUS_2000_I524911/d02-x01-y02
LogY=0
Title= Average $\cos 2\phi$ vrs $P_i$: 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0 GeV
XLabel=$P_i$ index
YLabel=$<\cos 2\phi>$
XMin=0
XMax=9.00
YMin=-0.10
YMax=0.15
YMajorTickMarks=10
YMinorTickMarks=3


