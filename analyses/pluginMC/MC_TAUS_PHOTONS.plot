BEGIN PLOT /MC_TAUS_PHOTONS/.*
LogY=1
END PLOT


BEGIN PLOT /MC_TAUS_PHOTONS/RestFramePhotonsEnergy.*
YLabel=$\frac{d\sigma}{dE} / \frac{\mathrm{pb}}{\mathrm{GeV}}$
XLabel=$\sum E_\gamma$ / GeV
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/RestFramePhotonsEnergyEl
Title=Photon energy in $\tau$ rest frame - electron decays
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/RestFramePhotonsEnergyMu
Title=Photon energy in $\tau$ rest frame - muonic decays
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/RestFramePhotonsEnergyHad
Title=Photon energy in $\tau$ rest frame - hadronic decays
END PLOT


BEGIN PLOT /MC_TAUS_PHOTONS/logPFracPhotons*
YLabel=$\frac{d\sigma}{d\log p_\mathrm{frac}} / \frac{\mathrm{pb}}{\mathrm{unit}}$
XLabel=$\log ( |\sum \vec p_\gamma| / |\vec p_\tau| )$
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/logPFracPhotonsEl
Title=Photon momentum fraction of $\tau$ decay - electronic decays
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/logPFracPhotonsMu
Title=Photon momentum fraction of $\tau$ decay - muonic decays
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/logPFracPhotonsHad
Title=Photon momentum fraction of $\tau$ decay - hadronic decays
END PLOT


BEGIN PLOT /MC_TAUS_PHOTONS/logPFracNotPhotons*
YLabel=$\frac{d\sigma}{d\log p_\mathrm{frac}} / \frac{\mathrm{pb}}{\mathrm{unit}}$
XLabel=$\log ( |\sum \vec p_{\neg \gamma}| / |\vec p_\tau| )$
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/logPFracPhotonsEl
Title=Non-photon momentum fraction of $\tau$ decay - electronic decays
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/logPFracPhotonsMu
Title=Non-hoton momentum fraction of $\tau$ decay - muonic decays
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/logPFracPhotonsHad
Title=Non-hoton momentum fraction of $\tau$ decay - hadronic decays
END PLOT


BEGIN PLOT /MC_TAUS_PHOTONS/nPhotons*
YLabel=$\frac{d\sigma}{dn_\gamma} / \frac{\mathrm{pb}}{\mathrm{unit}}$
XLabel=$n_\gamma$
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/nPhotonsEl
Title=Photon multiplicity in electronic $\tau$ decays
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/nPhotonsMu
Title=Photon multiplicity in muonic $\tau$ decays
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/nPhotonsHad
Title=Photon multiplicity in hadronic $\tau$ decays
END PLOT


BEGIN PLOT /MC_TAUS_PHOTONS/pFracPhotons*
YLabel=$\frac{d\sigma}{d\log p_\mathrm{frac}} / \frac{\mathrm{pb}}{\mathrm{unit}}$
XLabel=$|\sum \vec p_\gamma| / |\vec p_\tau|$
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/pFracPhotonsEl
Title=Photon momentum fraction of $\tau$ decay - electronic decays
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/pFracPhotonsMu
Title=Photon momentum fraction of $\tau$ decay - muonic decays
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/pFracPhotonsHad
Title=Photon momentum fraction of $\tau$ decay - hadronic decays
END PLOT


BEGIN PLOT /MC_TAUS_PHOTONS/tauMass*
YLabel=$\frac{d\sigma}{dm} / \frac{\mathrm{pb}}{GeV}$
XLabel=$m_\mathrm{decay\ products}$ / GeV
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/tauMassEl
Title=Invariant mass of $\tau$ decay products - electronic decays
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/tauMassMu
Title=Invariant mass of $\tau$ decay products - muonic decays
END PLOT

BEGIN PLOT /MC_TAUS_PHOTONS/tauMassHad
Title=Invariant mass of $\tau$ decay products - hadronic decays
END PLOT
