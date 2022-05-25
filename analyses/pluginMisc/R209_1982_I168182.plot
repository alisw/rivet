# BEGIN PLOT /R209_1982_I168182/DiMuon_pT
LegendYPos=0.8
LegendXPos=0.3
YMin=1E-5
YMax=10
Title=Drell-Yan production at $\sqrt{s}=62$ GeV
XLabel=$p_T$ (GeV)
YLabel=$1/\sigma d\sigma/dp_T$ (GeV$^{-1}$)
# + any additional plot settings you might like, see make-plots documentation
RatioPlot=0

# END PLOT
BEGIN PLOT /R209_1982_I168182/d0*
Title=Drell-Yan $\sqrt{s}=62$ GeV
LegendYPos=0.95
LegendXPos=0.2
END PLOT
BEGIN PLOT /R209_1982_I168182/d01-x01-y01
YMin=5E-5
YMax=5
XMax=20
XLabel=$M$ (GeV)
YLabel=$d\sigma/dM$ (bb/GeV)
# + any additional plot settings you might like, see make-plots documentation
END PLOT
BEGIN PLOT /R209_1982_I168182/d01-x01-y02
CustomLegend=$3.5 < M_{\mu^+\mu^-} < 18$ GeV
XLabel=$M$ (GeV)
YLabel=$d\sigma/dM$(nb/GeV)
# + any additional plot settings you might like, see make-plots documentation
END PLOT
BEGIN PLOT /R209_1982_I168182/d02-x01-y01
CustomLegend=$5 < M_{\mu^+\mu^-} < 8$ GeV
YMin=1E-4
YMax=5
XLabel=$p_T$ (GeV)
YLabel=$1/(2 p_T) d\sigma/dp_T$ (nb/GeV$^2$)
YMax=1
YMin=1E-4
#GofLegend=1
GofType=chi2
# + any additional plot settings you might like, see make-plots documentation
END PLOT

# ... add more histograms as you need them ...
