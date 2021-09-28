BEGIN PLOT /ARGUS_1992_I319102/d01-x01-y01
Title=$R=\sigma(e^+e^-\to \text{hadrons})/\sigma(e^+e^-\to \mu^+\mu^-)$
XLabel=$\sqrt{s}$/GeV
YLabel=$R$
LogY=0
ConnectGaps=1
END PLOT
BEGIN PLOT /ARGUS_1992_I319102/d02-x01-y01
Title=Charged Multiplicity at 10.47 GeV
XLabel=$N$
YLabel=Probability/\%
LogY=0
ConnectGaps=1
END PLOT
BEGIN PLOT /ARGUS_1992_I319102/d03-x01-y01
Title=Charged Multiplicity in $\Upsilon(4S)$ decays
XLabel=$N$
YLabel=Probability/\%
LogY=0
ConnectGaps=1
END PLOT
BEGIN PLOT /ARGUS_1992_I319102/sigma_hadrons
Title=$\sigma(e^+e^-\to \text{hadrons})$
XLabel=$\sqrt{s}$/GeV
YLabel=$\sigma(e^+e^-\to \text{hadrons})/pb$
LogY=0
ConnectGaps=1
XMin=9.
XMax=10.
END PLOT
BEGIN PLOT /ARGUS_1992_I319102/sigma_muons
Title=$\sigma(e^+e^-\to \mu^+\mu^-)$
XLabel=$\sqrt{s}$/GeV
YLabel=$\sigma(e^+e^-\to \mu^+\mu^-)/pb$
LogY=0
ConnectGaps=1
XMin=9.
XMax=10.
END PLOT

BEGIN PLOT /ARGUS_1992_I319102/d04-x01-y01
Title=Average charged Multiplicity at 10.47 GeV
XLabel=$\sqrt{s}$ [GeV]
YLabel=$\langle N_{\text{charged}}\rangle$
LogY=0
LegendYPos=0.2
END PLOT
BEGIN PLOT /ARGUS_1992_I319102/d05-x01-y01
Title=Average charged Multiplicity in $\Upsilon(4S)$ decays
XLabel=$\sqrt{s}$ [GeV]
YLabel=$\langle N_{\text{charged}}\rangle$
LogY=0
LegendYPos=0.2
END PLOT
