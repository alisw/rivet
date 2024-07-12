# BEGIN PLOT /MC_HHJETS/H_jet1_dR
Title=Separation between Higgs boson and leading jet
XLabel=$\Delta{R}(H, \mathrm{1st~jet})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta{R}(H, \mathrm{1st~jet})$ [pb]
RatioPlotYMin=0.1
RatioPlotYMax=1.8
# END PLOT

# BEGIN PLOT /MC_HHJETS/HH_dR
Title=Separation between Higgs bosons
XLabel=$\Delta{R}(H, H)$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta{R}(H, H)$ [pb]
#Rebin=2
RatioPlotYMin=0.1
RatioPlotYMax=2.0
# END PLOT

# BEGIN PLOT /MC_HHJETS/HH_dPhi
Title=Separation between Higgs bosons in $\phi$
XLabel=$\Delta\phi(H, H)$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta\phi(H, H)$ [pb]
#Rebin=2
RatioPlotYMin=0.1
RatioPlotYMax=2.0
# END PLOT

# BEGIN PLOT /MC_HHJETS/HH_deta
Title=Separation between Higgs bosons in $\eta$
XLabel=$\Delta\eta(H, H)$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta\eta(H, H)$ [pb]
#Rebin=2
RatioPlotYMin=0.1
RatioPlotYMax=2.0
# END PLOT

# BEGIN PLOT /MC_HHJETS/HH_pT$
Title=Di-Higgs $p_\perp$
XLabel=$p_\perp^{HH}$ [$\GeV$]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_\perp^{HH}$ [pb/$\GeV$]
LogX=0
#Rebin=4
#XMax=900
RatioPlotYMin=0.1
RatioPlotYMax=1.8
# END PLOT


# BEGIN PLOT /MC_HHJETS/H_jet1_deta
Title=Separation in $\eta$ between Higgs boson and leading jet
XLabel=$\Delta{\eta}(H, \mathrm{1st~jet})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta{\eta}(H, \mathrm{1st~jet})$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/HH_mass
Title=di-Higgs mass
XLabel=$m_{HH}$ [$\GeV$]
YLabel=$\mathrm{d}\sigma/\mathrm{d}m_{HH}$ [pb/$\GeV$]
#Rebin=2
#XMax=2E+03
# END PLOT

# BEGIN PLOT /MC_HHJETS/H_pT$
Title=Higgs boson $p_\perp$ (any Higgs)
XLabel=$p_\perp^H$ [$\GeV$]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_\perp^H$ [pb/$\GeV$]
LogX=0
#Rebin=4
#XMax=900
LegendXPos=0.15
LegendYPos=0.5
# END PLOT

# BEGIN PLOT /MC_HHJETS/H_pT1$
Title=Higgs boson $p_\perp$ (hardest Higgs)
XLabel=$p_\perp^H$ [$\GeV$]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_\perp^H$ [pb/$\GeV$]
LogX=0
#Rebin=4
#XMax=900
LegendXPos=0.15
LegendYPos=0.5
# END PLOT

# BEGIN PLOT /MC_HHJETS/H_pT2$
Title=Higgs boson $p_\perp$ (second hardest Higgs)
XLabel=$p_\perp^H$ [$\GeV$]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_\perp^H$ [pb/$\GeV$]
LogX=0
#Rebin=4
#XMax=900
LegendXPos=0.15
LegendYPos=0.5
# END PLOT


# BEGIN PLOT /MC_HHJETS/H_eta
Title=Higgs boson pseudorapidity
XLabel=$\eta_H$
YLabel=$\mathrm{d}\sigma/\mathrm{d} \eta_H$ [pb]
LegendXPos=0.35
LegendYPos=0.5
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/H_eta1
Title=Higgs boson pseudorapidity (hardest)
XLabel=$\eta_H$
YLabel=$\mathrm{d}\sigma/\mathrm{d} \eta_H$ [pb]
LegendXPos=0.35
LegendYPos=0.5
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/H_eta2
Title=Higgs boson pseudorapidity (second hardest)
XLabel=$\eta_H$
YLabel=$\mathrm{d}\sigma/\mathrm{d} \eta_H$ [pb]
LegendXPos=0.35
LegendYPos=0.5
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/H_phi
Title=Higgs boson azimuthal angle
XLabel=$\phi_H$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\phi_H$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/lepton_pT
Title=Lepton $p_\perp$
XLabel=$p_\perp^\ell$ [$\GeV$]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_\perp^\ell$ [pb/$\GeV$]
LogX=1
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/lepton_eta
Title=Lepton $\eta$
XLabel=$\eta_\ell$ [$\GeV$]
YLabel=$\mathrm{d}\sigma/\mathrm{d}\eta_\ell$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jets_dR_
Title=$\Delta R$ separation between jets
LegendXPos=0.10
LegendYPos=0.5
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jets_dR_12
XLabel=$\Delta{R}(\mathrm{jet~1}, \mathrm{jet~2})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta{R}(\mathrm{jet~1}, \mathrm{jet~2})$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jets_dR_13
XLabel=$\Delta{R}(\mathrm{jet~1}, \mathrm{jet~3})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta{R}(\mathrm{jet~1}, \mathrm{jet~3})$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jets_dR_23
XLabel=$\Delta{R}(\mathrm{jet~2}, \mathrm{jet~3})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta{R}(\mathrm{jet~2}, \mathrm{jet~3})$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jets_deta_
Title=Pseudorapidity separation between jets
LegendYPos=0.5
LegendXPos=0.30
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jets_deta_12
XLabel=$\Delta\eta(\mathrm{jet~1}, \mathrm{jet~2})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta\eta(\mathrm{jet~1}, \mathrm{jet~2})$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jets_deta_13
XLabel=$\Delta\eta(\mathrm{jet~1}, \mathrm{jet~3})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta\eta(\mathrm{jet~1}, \mathrm{jet~3})$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jets_deta_23
XLabel=$\Delta\eta(\mathrm{jet~2}, \mathrm{jet~3})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta\eta(\mathrm{jet~2}, \mathrm{jet~3})$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jets_dphi_
Title=Pseudorapidity separation between jets
LegendXPos=0.1
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jets_dphi_12
XLabel=$\Delta\phi(\mathrm{jet~1}, \mathrm{jet~2})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta\phi(\mathrm{jet~1}, \mathrm{jet~2})$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jets_dphi_13
XLabel=$\Delta\phi(\mathrm{jet~1}, \mathrm{jet~3})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta\phi(\mathrm{jet~1}, \mathrm{jet~3})$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jets_dphi_23
XLabel=$\Delta\phi(\mathrm{jet~2}, \mathrm{jet~3})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta\phi(\mathrm{jet~2}, \mathrm{jet~3})$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_eta_1
Title=Pseudorapidity of leading jet
XLabel=$\eta(\mathrm{jet~1})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\eta(\mathrm{jet~1})$ [pb]
LegendYPos=0.5
LegendXPos=0.30
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_eta_2
Title=Pseudorapidity of second jet
XLabel=$\eta(\mathrm{jet~2})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\eta(\mathrm{jet~2})$ [pb]
LegendYPos=0.5
LegendXPos=0.30
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_eta_3
Title=Pseudorapidity of third jet
XLabel=$\eta(\mathrm{jet~3})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\eta(\mathrm{jet~3})$ [pb]
LegendYPos=0.5
LegendXPos=0.30
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_eta_4
Title=Pseudorapidity of fourth jet
XLabel=$\eta(\mathrm{jet~4})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\eta(\mathrm{jet~4})$ [pb]
LegendYPos=0.5
LegendXPos=0.30
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_eta_pmratio_1
Title=Pseudorapidity $+/-$ ratio of first jet
XLabel=$\eta(\mathrm{jet~1})_+/\eta(\mathrm{jet~1})_-$
YLabel=$|\eta(\mathrm{jet~1}|$
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_eta_pmratio_2
Title=Pseudorapidity $+/-$ ratio of second jet
XLabel=$\eta(\mathrm{jet~2})_+/\eta(\mathrm{jet~2})_-$
YLabel=$|\eta(\mathrm{jet~2}|$
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_eta_pmratio_3
Title=Pseudorapidity $+/-$ ratio of third jet
XLabel=$\eta(\mathrm{jet~3})_+/\eta(\mathrm{jet~3})_-$
YLabel=$|\eta(\mathrm{jet~3}|$
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_eta_pmratio_4
Title=Pseudorapidity $+/-$ ratio of fourth jet
XLabel=$\eta(\mathrm{jet~4})_+/\eta(\mathrm{jet~4})_-$
YLabel=$|\eta(\mathrm{jet~4}|$
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_y_1
Title=Rapidity of first jet
XLabel=$y(\mathrm{jet~1})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}y(\mathrm{jet~1})$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_y_2
Title=Rapidity of second jet
XLabel=$y(\mathrm{jet~2})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}y(\mathrm{jet~2})$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_y_3
Title=Rapidity of third jet
XLabel=$y(\mathrm{jet~3})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}y(\mathrm{jet~3})$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_y_4
Title=Rapidity of fourth jet
XLabel=$y(\mathrm{jet~4})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}y(\mathrm{jet~4})$ [pb]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_y_pmratio_1
Title=Rapidity $+/-$ ratio of first jet
XLabel=$y(\mathrm{jet~1})_+/y(\mathrm{jet~1})_-$
YLabel=$|y(\mathrm{jet~1}|$
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_y_pmratio_2
Title=Rapidity $+/-$ ratio of second jet
XLabel=$y(\mathrm{jet~2})_+/y(\mathrm{jet~2})_-$
YLabel=$|y(\mathrm{jet~2}|$
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_y_pmratio_3
Title=Rapidity $+/-$ ratio of third jet
XLabel=$y(\mathrm{jet~3})_+/y(\mathrm{jet~3})_-$
YLabel=$|y(\mathrm{jet~3}|$
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_y_pmratio_4
Title=Rapidity $+/-$ ratio of fourth jet
XLabel=$y(\mathrm{jet~4})_+/y(\mathrm{jet~4})_-$
YLabel=$|y(\mathrm{jet~4}|$
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_mass_1
Title=Mass of first jet
XLabel=$m(\mathrm{jet~1})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}m(\mathrm{jet~1})$ [pb/$\GeV$]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_mass_2
Title=Mass of second jet
XLabel=$m(\mathrm{jet~2})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}m(\mathrm{jet~2})$ [pb/$\GeV$]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_mass_3
Title=Mass of third jet
XLabel=$m(\mathrm{jet~3})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}m(\mathrm{jet~3})$ [pb/$\GeV$]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_mass_4
Title=Mass of fourth jet
XLabel=$m(\mathrm{jet~4})$
YLabel=$\mathrm{d}\sigma/\mathrm{d}m(\mathrm{jet~4})$ [pb/$\GeV$]
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_multi_exclusive
Title=Exclusive jet multiplicity
XLabel=$N_{\mathrm{jet}}$
YLabel=$\sigma(N_{\mathrm{jet}})$ [pb]
XMajorTickMarks=10
XMinorTickMarks=0
ErrorBands=1
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_multi_inclusive
Title=Inclusive jet multiplicity
XLabel=$N_{\mathrm{jet}}$
YLabel=$\sigma(\geq N_{\mathrm{jet}})$ [pb]
XMajorTickMarks=10
XMinorTickMarks=0
ErrorBands=1
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_multi_ratio
Title=Ratio of jet multiplicity
XLabel=$N_{\mathrm{jet}}$
YLabel=$\sigma(\geq N_{\mathrm{jet}})/\sigma(\geq N_{\mathrm{jet}}-1)$
XMajorTickMarks=10
XMinorTickMarks=0
LogY=0
ErrorBands=1
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/log10_R_0
Title=$\log_{10}$(Integrated $0$ jet rate in $k_\perp$ [$\GeV$])
XLabel=$\log_{10}(d_{\mathrm{cut}}/\GeV)$
YLabel=$R_{0}$
#Rebin=2
LegendYPos=0.8
LegendXPos=0.75
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/log10_R_1
Title=$\log_{10}$(Integrated $1$ jet rate in $k_\perp$ [$\GeV$])
XLabel=$\log_{10}(d_{\mathrm{cut}}/\GeV)$
YLabel=$R_{1}$
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/log10_R_2
Title=$\log_{10}$(Integrated $2$ jet rate in $k_\perp$ [$\GeV$])
XLabel=$\log_{10}(d_{\mathrm{cut}}/\GeV)$
YLabel=$R_{2}$
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/log10_R_3
Title=$\log_{10}$(Integrated $3$ jet rate in $k_\perp$ [$\GeV$])
XLabel=$\log_{10}(d_{\mathrm{cut}}/\GeV)$
YLabel=$R_{3}$
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/log10_R_4
Title=$\log_{10}$(Integrated $4$ jet rate in $k_\perp$ [$\GeV$])
XLabel=$\log_{10}(d_{\mathrm{cut}}/\GeV)$
YLabel=$R_{\geq4}$
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/log10_d_01
Title=$\log_{10}$($k_\perp$ jet resolution $0 \to 1$ [$\GeV$])
XLabel=$\log_{10}(d_{01}/\GeV)$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\log_{10}(d_{01}/\GeV)$ [pb]
LegendXPos=0.15
LegendYPos=0.5
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/log10_d_12
Title=$\log_{10}$($k_\perp$ jet resolution $1 \to 2$ [$\GeV$])
XLabel=$\log_{10}(d_{12}/\GeV)$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\log_{10}(d_{12}/\GeV)$ [pb]
LegendXPos=0.15
LegendYPos=0.5
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/log10_d_23
Title=$\log_{10}$($k_\perp$ jet resolution $2 \to 3$ [$\GeV$])
XLabel=$\log_{10}(d_{23}/\GeV)$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\log_{10}(d_{23}/\GeV)$ [pb]
LegendXPos=0.15
LegendYPos=0.5
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/log10_d_34
Title=$\log_{10}$($k_\perp$ jet resolution $3 \to 4$ [$\GeV$])
XLabel=$\log_{10}(d_{34}/\GeV)$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\log_{10}(d_{34}/\GeV)$ [pb]
LegendXPos=0.15
LegendYPos=0.5
#Rebin=2
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_pT_1
Title=Transverse momentum of leading jet
XLabel=$p_\perp(\mathrm{jet~1})$ [$\GeV$]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_\perp(\mathrm{jet~1})$ [pb/$\GeV$]
LogX=1
LegendXPos=0.05
LegendYPos=0.5
#XMin=20
#XMax=1500.0
RatioPlotYMin=0.1
RatioPlotYMax=1.6
#Rebin=8
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_pT_2
Title=Transverse momentum of second jet
XLabel=$p_\perp(\mathrm{jet~2})$ [$\GeV$]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_\perp(\mathrm{jet~2})$ [pb/$\GeV$]
LogX=1
LegendXPos=0.05
LegendYPos=0.5
#XMin=20.0
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_pT_3
Title=Transverse momentum of third jet
XLabel=$p_\perp(\mathrm{jet~3})$ [$\GeV$]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_\perp(\mathrm{jet~3})$ [pb/$\GeV$]
LogX=1
LegendXPos=0.05
LegendYPos=0.5
#XMin=20.0
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_pT_4
Title=Transverse momentum of fourth jet
XLabel=$p_\perp(\mathrm{jet~4})$ [$\GeV$]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_\perp(\mathrm{jet~4})$ [pb/$\GeV$]
LogX=1
LegendXPos=0.05
LegendYPos=0.5
#XMin=20.0
# END PLOT

# BEGIN PLOT /MC_HHJETS/jet_HT
Title=Scalar sum of jet transverse momenta ($H_T$)
XLabel=$H_T$ [$\GeV$]
YLabel=$\mathrm{d}\sigma/\mathrm{d}H_T$ [pb/$\GeV$]
LogX=1
# END PLOT

# BEGIN PLOT /MC_HHJETS/jets_mjj
Title=Dijet invariant mass spectrum
XLabel=$m_{jj}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}m_{jj}$ [pb/GeV]
LegendAlign=r
# END PLOT
