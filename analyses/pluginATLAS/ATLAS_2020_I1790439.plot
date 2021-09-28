BEGIN PLOT /ATLAS_2020_I1790439/d..
LegendAlign=r
XTwosidedTicks=1
YTwosidedTicks=1
LogY=0
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d03-x01-y01
LegendAlign=l
LegendXPos=0.05
YLabel=Cross section [fb]
XCustomMajorTicks=1  $\sigma_{4\mu}$  2  $\sigma_{4e}$  3  $\sigma_{2\mu2e}$  4  $\sigma_{2e2\mu}$  5  $\sigma_{4\mu+4e}$  6  $\sigma_{2\cdot(2\mu2e)}$  7  $\sigma_\text{sum}$  8  $\sigma_\text{comb}$  9  $\sigma_\text{tot}$
XMinorTickMarks=0
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d05-x01-y01
XLabel=$p_\text{T}^{4l}$ [GeV]
YLabel=d$\sigma$/d$p_\text{T}^{4l}$ [fb/GeV]
XMin=1
LogX=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d07-x01-y01
LegendAlign=l
LegendXPos=0.05
XLabel=$m_{12}$ [GeV]
YLabel=d$\sigma$/d$m_{12}$ [fb/GeV]
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d09-x01-y01
XLabel=$m_{34}$ [GeV]
YLabel=d$\sigma$/d$m_{34}$ [fb/GeV]
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d11-x01-y01
XLabel=$|y_{4\ell}|$
YLabel=d$\sigma$/d$|y_{4\ell}|$ [fb]
RatioPlotYMax=1.99
RatioPlotYMin=0.01
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d13-x01-y01
XLabel=$|\cos(\theta^*)|$
YLabel=d$\sigma$/d$|\cos(\theta^*)|$ [fb]
RatioPlotYMax=2.19
RatioPlotYMin=0.0
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d15-x01-y01
XLabel=$\cos(\theta_1)$
YLabel=d$\sigma$/d$\cos(\theta_1)$ [fb]
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d17-x01-y01
XLabel=$\cos(\theta_2)$
YLabel=d$\sigma$/d$\cos(\theta_2)$ [fb]
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d19-x01-y01
XLabel=$\phi$
YLabel=d$\sigma$/d$\phi$ [fb]
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d21-x01-y01
XLabel=$\phi_1$
YLabel=d$\sigma$/d$\phi_1$ [fb]
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d23-x01-y01
XLabel=$N_\text{jets}$
YLabel=d$\sigma$/d$N_\text{jets}$ [fb]
XMinorTickMarks=0
XCustomMajorTicks=0.0  $=0$  1.0  $=1$  2.0  $=2$  3.0  $\geq3$
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d25-x01-y01
XLabel=$N_\text{jets}$
YLabel=d$\sigma$/d$N_\text{jets}$ [fb]
XMinorTickMarks=0
XCustomMajorTicks=0.0  $\geq0$  1.0  $\geq1$  2.0  $\geq2$  3.0  $\geq3$
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d26-x01-y01
YLabel=d$\sigma$/d$N_{b-\text{jets}}$ [fb]
XMinorTickMarks=0
XCustomMajorTicks=1.0  $N_\text{jets}=0$  2.0  $N_{b-\text{jets}}=0$  3.0  $N_{b-\text{jets}}\geq1$
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d28-x01-y01
XLabel=leading jet $p_\text{T}$ [GeV]
YLabel=d$\sigma$/d$p_\text{T}$ [fb/GeV]
XCustomMajorTicks=15.  $N_\text{jets}=0\qquad\qquad$  30.  30  60.  60  120.  120  350.  350
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d30-x01-y01
XLabel=subleading jet $p_\text{T}$ [GeV]
YLabel=d$\sigma$/d$p_\text{T}$ [fb/GeV]
XCustomMajorTicks=15.  $N_\text{jets}\leq1\qquad\qquad$  30.  30  60.  60  120.  120  350.  350
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d32-x01-y01
XLabel=$m_{jj}$ [GeV]
YLabel=d$\sigma$/d$m_{jj}$ [fb/GeV]
XCustomMajorTicks=-250.  $N_\text{jets}=0\qquad\qquad$  60.  $60\quad$  120.  $\quad120$  450.  450  3000.  3000
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d34-x01-y01
XLabel=$\Delta\eta_{jj}$
YLabel=d$\sigma$/d$\Delta\eta_{jj}$ [fb]
XCustomMajorTicks=-0.5  $N_\text{jets}\leq1\qquad\qquad$  0.  0  1.  1  2.5  2.5  9.  9
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d36-x01-y01
XLabel=$\Delta\phi_{jj}$
YLabel=d$\sigma$/d$\Delta\phi_{jj}$ [fb]
XCustomMajorTicks=-0.5  $N_\text{jets}\leq1\qquad\qquad$  0.  0  1.57  $\pi/2$  3.14  $\pi$  4.71  $3\pi/2$  6.28  $2\pi$
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d38-x01-y01
XLabel=$p_\text{T}^{4\ell j}$ [GeV]
YLabel=d$\sigma$/d$p_\text{T}^{4\ell j}$ [fb/GeV]
XCustomMajorTicks=-30.  $N_\text{jets}=0\qquad\qquad$  0.  0  60.  60  120  120  350  350
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d40-x01-y01
XLabel=$p_\text{T}^{4\ell jj}$ [GeV]
YLabel=d$\sigma$/d$p_\text{T}^{4\ell jj}$ [fb/GeV]
XCustomMajorTicks=-30.  $N_\text{jets}\leq1\qquad\qquad$  0.  0  60.  60  120  120  350  350
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d42-x01-y01
XLabel=$m_{4\ell j}$ [GeV]
YLabel=d$\sigma$/d$m_{4\ell j}$ [fb/GeV]
XCustomMajorTicks=60.  $N_\text{jets}=0\qquad\qquad$  120  120  250  250  400  400  600  600  2000  2000
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d44-x01-y01
XLabel=$m_{4\ell jj}$ [GeV]
YLabel=d$\sigma$/d$m_{4\ell jj}$ [fb/GeV]
XCustomMajorTicks=60.  $N_\text{jets}\leq1\qquad\qquad$  180  180  400  400  600  600  1000  1000  2500  2500
LogY=1
RatioPlotYMax=1.7
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d46-x01-y01
Title=$H\to 4\ell$
XLabel=$m_{12}$ vs $m_{34}$
YLabel=d$\sigma$/d$m_{12}$d$m_{34}$ [fb]
XMinorTickMarks=0
LegendAlign=l
LegendXPos=0.05
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d48-x01-y01
Title=$H\to 2\ell 2\mu$
XLabel=$m_{12}$ vs $m_{34}$
YLabel=d$\sigma$/d$m_{12}$d$m_{34}$ [fb]
XMinorTickMarks=0
LegendAlign=l
LegendXPos=0.05
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d49-x01-y01
Title=$H\to 2\ell 2e$
XLabel=$m_{12}$ vs $m_{34}$
YLabel=d$\sigma$/d$m_{12}$d$m_{34}$ [fb]
XMinorTickMarks=0
LegendAlign=l
LegendXPos=0.05
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d51-x01-y01
Title=$0.0<|y_{4\ell}|<0.5$
XLabel=$p_\text{T}^{4\ell}$
YLabel=d$\sigma$/d$p_{T}^{4\ell}$ [fb]
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d51-x01-y03
Title=$0.5<|y_{4\ell}|<1.0$
XLabel=$p_\text{T}^{4\ell}$
YLabel=d$\sigma$/d$p_{T}^{4\ell}$ [fb]
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d51-x01-y05
Title=$1.0<|y_{4\ell}|<1.5$
XLabel=$p_\text{T}^{4\ell}$
YLabel=d$\sigma$/d$p_{T}^{4\ell}$ [fb]
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d51-x01-y07
Title=$1.5<|y_{4\ell}|<2.5$
XLabel=$p_\text{T}^{4\ell}$
YLabel=d$\sigma$/d$p_{T}^{4\ell}$ [fb]
LogY=1
RatioPlotYMax=1.9
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d53-x01-y01
Title=$N_\text{jets}=0$
XLabel=$p_\text{T}^{4\ell}$
YLabel=d$\sigma$/d$p_{T}^{4\ell}$ [fb]
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d53-x01-y03
Title=$N_\text{jets}=1$
XLabel=$p_\text{T}^{4\ell}$
YLabel=d$\sigma$/d$p_{T}^{4\ell}$ [fb]
LogY=1
RatioPlotYMax=1.9
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d53-x01-y05
Title=$N_\text{jets}=2$
XLabel=$p_\text{T}^{4\ell}$
YLabel=d$\sigma$/d$p_{T}^{4\ell}$ [fb]
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d53-x01-y07
Title=$N_\text{jets}\geq 3$
XLabel=$p_\text{T}^{4\ell}$
YLabel=d$\sigma$/d$p_{T}^{4\ell}$ [fb]
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d55-x01-y01
XLabel=$p_\text{T}^{4\ell}$ vs $p_\text{T}^{4\ell j}$
YLabel=d$^2\sigma$/d$p_\text{T}^{4\ell}$d$p_\text{T}^{4\ell j}$ [fb]
XMinorTickMarks=0
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d57-x01-y01
XLabel=$p_\text{T}^{4\ell}$ vs $m_{4\ell j}$
YLabel=d$^2\sigma$/d$p_\text{T}^{4\ell}$d$m_{4\ell j}$ [fb]
XMinorTickMarks=0
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d59-x01-y01
XLabel=$p_\text{T}^{4\ell}$ vs $p_\text{T}^\text{jet1}$
YLabel=d$^2\sigma$/d$p_\text{T}^{4\ell}$d$p_\text{T}^\text{jet1}$ [fb]
XMinorTickMarks=0
RatioPlotYMin=0.3
RatioPlotYMax=1.7
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d61-x01-y01
XLabel=$p_\text{T}^{\text{jet1}}$ vs $|y^\text{jet1}|$
YLabel=d$^2\sigma$/d$p_\text{T}^{\text{jet1}}$d$|y^\text{jet1}|$ [fb]
XMinorTickMarks=0
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d63-x01-y01
XLabel=$p_\text{T}^{\text{jet1}}$ vs $p_\text{T}^\text{jet2}$
YLabel=d$^2sigma$/d$p_\text{T}^{\text{jet1}}$d$p_\text{T}^\text{jet2}$ [fb]
XMinorTickMarks=0
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d65-x01-y01
Title=$H\to 4\ell$
XLabel=$m_{12}$ [GeV]
YLabel=d$\sigma$/d$m_{12}$ [fb/GeV]
LegendAlign=l
LegendXPos=0.05
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d66-x01-y01
Title=$H\to 2\ell 2\ell^\prime$
XLabel=$m_{12}$ [GeV]
YLabel=d$\sigma$/d$m_{12}$ [fb/GeV]
LegendAlign=l
LegendXPos=0.05
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d68-x01-y01
Title=$H\to 4\ell$
XLabel=$m_{34}$ [GeV]
YLabel=d$\sigma$/d$m_{34}$ [fb/GeV]
LogY=1
RatioPlotYMax=1.7
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d69-x01-y01
Title=$H\to 2\ell 2\ell^\prime$
XLabel=$m_{34}$ [GeV]
YLabel=d$\sigma$/d$m_{34}$ [fb/GeV]
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d71-x01-y01
Title=$H\to 4\ell$
XLabel=$\phi$
YLabel=d$\sigma$/d$\phi$ [fb]
LogY=1
RatioPlotYMax=1.9
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d72-x01-y01
Title=$H\to 2\ell 2\ell^\prime$
XLabel=$\phi$
YLabel=d$\sigma$/d$\phi$ [fb]
LogY=1
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d74-x01-y01
Title=$H\to 4\ell$
XLabel=$m_{12}$ vs $m_{34}$
YLabel=d$^2\sigma$/d$m_{12}$d$m_{34}$ [fb]
LogY=1
XMinorTickMarks=0
LegendAlign=l
LegendXPos=0.05
END PLOT

BEGIN PLOT /ATLAS_2020_I1790439/d75-x01-y01
Title=$H\to 2\ell 2\ell^\prime$
XLabel=$m_{12}$ vs $m_{34}$
YLabel=d$^2\sigma$/d$m_{12}$d$m_{34}$ [fb]
LogY=1
XMinorTickMarks=0
LegendAlign=l
LegendXPos=0.05
END PLOT

