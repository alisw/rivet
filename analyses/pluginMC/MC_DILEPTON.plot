BEGIN PLOT /MC_DILEPTON/com_costheta.*
Title=Angle between leptons and boost in boson frame
XLabel=$\cos\theta^*$
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}\cos\theta^*$
LegendYPos=0.4
END PLOT

BEGIN PLOT /MC_DILEPTON/com_costheta_l1
Title=Angle between leading lepton and boost in boson frame
END PLOT

BEGIN PLOT /MC_DILEPTON/com_costheta_l2
Title=Angle between subleading lepton and boost in boson frame
LegendYPos=0.9
END PLOT


BEGIN PLOT /MC_DILEPTON/com_ppara.*
Title=Lepton momenta along boost in boson frame
XLabel=$p^*_\parallel$ [GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}p^*_\parallel$ [1/GeV]
LegendYPos=0.4
END PLOT

BEGIN PLOT /MC_DILEPTON/com_ppara_l1
Title=Leading lepton momentum along boost in boson frame
END PLOT

BEGIN PLOT /MC_DILEPTON/com_ppara_l2
Title=Subleading lepton momentum along boost in boson frame
LegendYPos=0.9
END PLOT


BEGIN PLOT /MC_DILEPTON/com_pperp
Title=Lepton momenta transverse to boost in boson frame
XLabel=$p^*_\perp$ [GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}p^*_\perp$ [1/GeV]
END PLOT


BEGIN PLOT /MC_DILEPTON/lep._pt
XLabel=$p_\mathrm{T}$ [GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}p_\mathrm{T}$ [1/GeV]
END PLOT

BEGIN PLOT /MC_DILEPTON/lep1_pt
Title=Leading lepton transverse momentum in lab frame
END PLOT

BEGIN PLOT /MC_DILEPTON/lep2_pt
Title=Subleading lepton transverse momentum in lab frame
END PLOT


BEGIN PLOT /MC_DILEPTON/lep._costheta
XLabel=$\cos\theta$
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}\cos\theta$
LegendXPos=0.4
END PLOT

BEGIN PLOT /MC_DILEPTON/lep1_costheta
Title=Angle between leading lepton and boost in lab frame
END PLOT

BEGIN PLOT /MC_DILEPTON/lep2_costheta
Title=Angle between subleading lepton and boost in lab frame
END PLOT


BEGIN PLOT /MC_DILEPTON/lep._ppara
XLabel=$p_\parallel$ [GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}p_\parallel$ [1/GeV]
LegendYPos=0.4
END PLOT

BEGIN PLOT /MC_DILEPTON/lep1_ppara
Title=Leading lepton momentum along boost in lab frame
END PLOT

BEGIN PLOT /MC_DILEPTON/lep2_ppara
Title=Subleading lepton momentum along boost in lab frame
END PLOT


BEGIN PLOT /MC_DILEPTON/lep._pperp
XLabel=$p_\perp$ [GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}p_\perp$ [1/GeV]
END PLOT

BEGIN PLOT /MC_DILEPTON/lep1_pperp
Title=Leading lepton momentum transverse to boost in lab frame
END PLOT

BEGIN PLOT /MC_DILEPTON/lep2_pperp
Title=Subleading lepton momentum transverse to boost in lab frame
END PLOT
