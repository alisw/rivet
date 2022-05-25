BEGIN PLOT /H1_1995_I394793/d01-x01-y01
Title= Cos theta with low Q
XLabel= $\cos \Theta_B$
YLabel= $1/N dn_{ch} / d \cos \Theta_B $
# + any additional plot settings you might like, see make-plots documentation

END PLOT

# ... add more histograms as you need them ...
BEGIN PLOT /H1_1995_I394793/d01-x01-y02
Title= Cos theta with high Q
XLabel= $\cos \Theta_B$
YLabel= $1/N dn_{ch} / d \cos \Theta_B $
END PLOT

BEGIN PLOT /H1_1995_I394793/d02-x01-y01
Title= Cos theta with low Q (no Efwd)
XLabel= $\cos \Theta_B$
YLabel= $1/N dn_{ch} / d \cos \Theta_B $
END PLOT

BEGIN PLOT /H1_1995_I394793/d02-x01-y02
Title= Cos theta with high Q (no Efwd)
XLabel= $\cos \Theta_B$
YLabel= $1/N dn_{ch} / d \cos \Theta_B $
END PLOT

BEGIN PLOT /H1_1995_I394793/d03-x01-y01
Title= $D( x_p)$ positive, low Q
XLabel= $x_p$
YLabel= $1/N  dn^+/dx_p$
END PLOT

BEGIN PLOT /H1_1995_I394793/d03-x01-y02
Title= $D( x_p)$ negative, low Q
XLabel= $x_p$
YLabel= $1/N  dn^-/dx_p$
END PLOT

BEGIN PLOT /H1_1995_I394793/d03-x01-y03
Title= $D( x_p)$ positive, high Q
XLabel= $x_p$
YLabel= $1/N  dn^+/dx_p$
END PLOT

BEGIN PLOT/H1_1995_I394793/d03-x01-y04
Title= $D( x_p)$ negative, high Q
XLabel= $x_p$
YLabel= $1/N  dn^-/dx_p$
END PLOT

BEGIN PLOT/H1_1995_I394793/d04-x01-y01
Title= D(xi) the low Q2
XLabel= $\xi$
YLabel= $D(\xi)$
LogY=1
END PLOT

BEGIN PLOT /H1_1995_I394793/d04-x01-y02
Title= D(xi) the high Q2
XLabel= $\xi$
YLabel= $D(\xi)$
LogY=1
END PLOT

BEGIN PLOT /H1_1995_I394793/d05-x01-y01
Title= Average charged particle multiplicity vrs. Q2
XLabel=$Q^2$ [GeV$^2$]
YLabel=$<n_{ch}>$
LogX=1
#LegendYPos=0.6
LegendXPos=0.1
END PLOT

BEGIN PLOT /H1_1995_I394793/d06-x01-y01
Title= Average charged particle multiplicity vrs. Q2
XLabel=$Q^2$ [GeV$^2$]
YLabel=$<n_{ch}>$
LogX=1
END PLOT
