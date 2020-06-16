# BEGIN PLOT /MC_PARTONICTOPS/t_.*_n
XLabel=$n$
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}n$
# END PLOT

# BEGIN PLOT /MC_PARTONICTOPS/t_.*_pT
XLabel=$p_\perp$~[GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}p_\perp$~[GeV$^{-1}$]
# END PLOT

# BEGIN PLOT /MC_PARTONICTOPS/t_.*_y
XLabel=$y$
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}y$
# END PLOT


# BEGIN PLOT /MC_PARTONICTOPS/t_all_n
Title=Top-quark per-event multiplicity
# END PLOT
# BEGIN PLOT /MC_PARTONICTOPS/t_all_pT
Title=Top-quark transverse momentum spectrum
# END PLOT
# BEGIN PLOT /MC_PARTONICTOPS/t_all_y
Title=Top-quark rapidity distribution
# END PLOT

# BEGIN PLOT /MC_PARTONICTOPS/t_all_pT_dfirstlast$
Title=First-to-last top-quark transverse momentum difference
XLabel=$\Delta p_\perp = p_\perp^\mathrm{last} - p_\perp^\mathrm{first}$~[GeV]
YLabel=$1/\sigma \, \mathrm{d}\sigma/\mathrm{d}\Delta p_\perp$~[GeV$^{-1}$]
# END PLOT

# BEGIN PLOT /MC_PARTONICTOPS/t_all_pT_dfirstlast_prof
Title=First-to-last top-quark $p_\perp$ difference vs $p_\perp^\mathrm{first}$
XLabel=$p_\perp^\mathrm{first}$~[GeV]
YLabel=$\langle \Delta p_\perp \rangle$~[GeV]
LogY=0
LegendYPos=0.2
# END PLOT

# BEGIN PLOT /MC_PARTONICTOPS/t_lep_n
Title=Leptonic top-quark per-event multiplicity
# END PLOT
# BEGIN PLOT /MC_PARTONICTOPS/t_lep_pT
Title=Leptonic top-quark transverse momentum spectrum
# END PLOT
# BEGIN PLOT /MC_PARTONICTOPS/t_lep_y
Title=Leptonic top-quark rapidity distribution
# END PLOT

# BEGIN PLOT /MC_PARTONICTOPS/t_had_n
Title=Hadronic top-quark per-event multiplicity
# END PLOT
# BEGIN PLOT /MC_PARTONICTOPS/t_had_pT
Title=Hadronic top-quark transverse momentum spectrum
# END PLOT
# BEGIN PLOT /MC_PARTONICTOPS/t_had_y
Title=Hadronic top-quark rapidity distribution
# END PLOT
