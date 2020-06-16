BEGIN PLOT /HRS_1987_I250823/
LogY=0
END PLOT
BEGIN PLOT /HRS_1987_I250823/d0[1,2,3]-x0[1,2,3]-y01
XLabel=$x_E$
YLabel=$\rho_{00}$
Title=$\rho_{00}$ for $D^*$ production, helicity-beam frame
END PLOT
BEGIN PLOT /HRS_1987_I250823/d0[1,2,3]-x0[1,2,3]-y02
XLabel=$x_E$
YLabel=$\rho_{1-1}$
Title=$\rho_{1-1}$ for $D^*$ production, helicity-beam frame
END PLOT

BEGIN PLOT /HRS_1987_I250823/d0[1,2,3]-x0[1,2,3]-y03
XLabel=$x_E$
YLabel=$\text{Re}\rho_{10}$
Title=$\Text{Re}\rho_{10}$ for $D^*$ production, helicity-beam frame
END PLOT

BEGIN PLOT /HRS_1987_I250823/d04-x01-y01
YLabel=$\rho_{00}$
Title=$\rho_{00}$ for $D^*$ production, helicity-beam frame, ($x_E>0.4$)
END PLOT
BEGIN PLOT /HRS_1987_I250823/d04-x01-y02
YLabel=$\rho_{1-1}$
Title=$\rho_{1-1}$ for $D^*$ production, helicity-beam frame, ($x_E>0.4$)
END PLOT
BEGIN PLOT /HRS_1987_I250823/d04-x01-y03
YLabel=$\text{Re}\rho_{10}$
Title=$\Text{Re}\rho_{10}$ for $D^*$ production, helicity-beam frame, ($x_E>0.4$)
END PLOT

BEGIN PLOT /HRS_1987_I250823/d05-x01-y0[1,2,3]
XCustomMajorTicks=1.  All	2.  $p_\perp<0.75$~GeV	3.  $p_\perp>0.75$~GeV
END PLOT

BEGIN PLOT /HRS_1987_I250823/d06-x01-y0[1,2,3]
XCustomMajorTicks=1.  $S\geq0.1$		2.  $S<0.1$
END PLOT

BEGIN PLOT /HRS_1987_I250823/d0[5,6]-x01-y01
YLabel=$\rho_{00}$
Title=$\rho_{00}$ for $D^*$ production, helicity-quark frame, ($x_E>0.4$)
END PLOT
BEGIN PLOT /HRS_1987_I250823/d0[5,6]-x01-y02
YLabel=$\rho_{1-1}$
Title=$\rho_{1-1}$ for $D^*$ production, helicity-quark frame, ($x_E>0.4$)
END PLOT
BEGIN PLOT /HRS_1987_I250823/d0[5,6]-x01-y03
YLabel=$\text{Re}\rho_{10}$
Title=$\Text{Re}\rho_{10}$ for $D^*$ production, helicity-quark frame, ($x_E>0.4$)
END PLOT
