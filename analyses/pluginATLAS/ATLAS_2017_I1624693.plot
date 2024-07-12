BEGIN PLOT /ATLAS_2017_I1624693/d01-x01-y01
Title=
XLabel=$Q$ [GeV]
YLabel=$\Delta (Q)$
LogY=0
YMin=-0.002
LegendXPos=0.98
LegendAlign=r
RatioPlotMode=deviation
END PLOT

BEGIN PLOT /ATLAS_2017_I1624693/d02-x01-y01
Title=
XLabel=$Q$ [GeV]
YLabel=$\Delta_\mathrm{3h} (Q)$
LogY=0
LegendXPos=0.98
LegendYPos=0.3
LegendAlign=r
RatioPlotMode=deviation
END PLOT

BEGIN PLOT /ATLAS_2017_I1624693/d03-x01-y01
Title=$m_\mathrm{3h} < 0.59$ GeV
XLabel=$X$
YLabel=$Y$
ZLabel=$N_\text{3h} / N_\text{ch}$
ZLabelSep=6.5
RightMargin=1.0
LogZ=1
ZMax=0.005
ZMin=0.00000011
END PLOT

BEGIN PLOT /ATLAS_2017_I1624693/d03-x01-y02
Title=$m_\mathrm{3h} < 0.59$ GeV
XLabel=$X$
YLabel=$Y$
ZLabel=(data - MC)/$\sigma$
ZLabelSep=6.5
LogZ=0
ZMax=20
ZMin=-5
END PLOT

