# BEGIN PLOT /MC_XS/XS
Title=Total generated cross section
XLabel=
YLabel=$\sigma_\text{total}$ [pb]
LogY=0
ShowZero=0
XCustomMajorTicks=0.	$\quad$
#XMajorTickMarks=20
XMinorTickMarks=0
LegendAlign=r
# END PLOT

# BEGIN PLOT /MC_XS/N
Title=Number of generated events
XLabel=
YLabel=$N$
LogY=0
XCustomMajorTicks=0.5	$\quad$
#XMajorTickMarks=20
XMinorTickMarks=0
LegendXPos=0.05
LegendYPos=0.15
# END PLOT

# BEGIN PLOT /MC_XS/pmXS
Title=Fraction of positive and negative weighted events
XLabel=$\text{sgn}(w)$
XCustomMajorTicks=-0.5	$w<0$	0.5	$w\geq0$
YLabel=$\text{d}\sigma/\text{d sgn}(w)$ [pb]
LogY=0
ShowZero=0
XMinorTickMarks=0
XMin=-1
XMax=1
LegendXPos=0.05
# END PLOT

# BEGIN PLOT /MC_XS/pmN
Title=Number of positive and negative weighted events
XLabel=$\text{sgn}(w)$
YLabel=$N_\pm$
LogY=0
ShowZero=0
XCustomMajorTicks=-0.5	$w<0$	0.5	$w\geq0$
XMinorTickMarks=0
XMin=-1
XMax=1
LegendXPos=0.05
# END PLOT
