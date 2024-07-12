# BEGIN PLOT /MC_WEIGHTS/.*weight.*
Title=Event weight distribution
XLabel=$w$
YLabel=$1/N \mathrm{d}N / \mathrm{d}w$
LegendAlign=r
# END PLOT


# BEGIN PLOT /MC_WEIGHTS/xsfraction_neg
LegendAlign=r
Title=Negative weight fraction
XLabel=
YLabel=$\sum_{w_i<0} |w_i| / \sum |w_i|$
LogY=0
ShowZero=0
XCustomMajorTicks=0.5	$f_{w<0}$
XMinorTickMarks=0
# END PLOT
