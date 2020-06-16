BEGIN PLOT /E605_1991_I302822/d*
Title=E605 $\sqrt{s}=38.8 $ GeV, $ -0.1 < x_F < 0.2 $
LegendYPos=0.9
LegendXPos=0.4
RatioPlotYMax=4.
#RatioPlotMode=datamc
END PLOT

BEGIN PLOT /E605_1991_I302822/d17-x01-y01
XLabel=$p_T$
YLabel=$d\sigma/dp_T$
CustomLegend=$7 < p_T < 8 $ GeV
# + any additional plot settings you might like, see make-plots documentation
END PLOT
BEGIN PLOT /E605_1991_I302822/d18-x01-y01
XLabel=$p_T$
YLabel=$d\sigma/dp_T$
CustomLegend=$8 < M_{\ell\ell} < 9 $ GeV
# + any additional plot settings you might like, see make-plots documentation
END PLOT
BEGIN PLOT /E605_1991_I302822/d19-x01-y01
XLabel=$p_T$
YLabel=$d\sigma/dp_T$
CustomLegend=$10.5 < M_{\ell\ell} < 11.5 $ GeV
# + any additional plot settings you might like, see make-plots documentation
END PLOT
BEGIN PLOT /E605_1991_I302822/d20-x01-y01
XLabel=$p_T$
YLabel=$d\sigma/dp_T$
CustomLegend=$11.5 < M_{\ell\ell} < 13.5 $ GeV
# + any additional plot settings you might like, see make-plots documentation
END PLOT
BEGIN PLOT /E605_1991_I302822/d21-x01-y01
XLabel=$p_T$
YLabel=$d\sigma/dp_T$
CustomLegend=$13.5 < M_{\ell\ell} < 18 $ GeV
# + any additional plot settings you might like, see make-plots documentation
END PLOT

# ... add more histograms as you need them ...
