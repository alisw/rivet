## Flags for building all plugins
AM_LDFLAGS = $(LDFLAGS) -module -avoid-version -L$(FASTJETLIBPATH)
LIBS = $(FASTJETCONFIGLIBADD)
lib_LTLIBRARIES =

noinst_LTLIBRARIES = libRivetAnalysisTools.la
libRivetAnalysisTools_la_SOURCES = \
    MC_ParticleAnalysis.cc \
    MC_JetAnalysis.cc \
    MC_JetSplittings.cc

## ANALYSIS CATEGORIES
##
## Unvalidated analyses: add new standard analyses here, and only
## move them into the collider-specific standard plugin libraries
## once they have been finished and checked. The --enable-unvalidated
## flag is needed to install the unvalidated analyses.
##
## Preliminary analyses: validated analyses whose experimental paper has not
## been fully accepted for publication should go here. Analyses in this group
## are considered safe to use but the reference data may yet change. In
## progressing from preliminary status to a permanent published analyses
## collection, the analysis name is very likely to change, so you should ensure
## that any Rivet-using scripts are not broken by such name changes when
## upgrading between Rivet versions. These analyses will not be available if
## Rivet is built with the --disable-preliminary configure flag.
##
## Obsolete analyses: as mentioned above, when a preliminary analysis becomes
## permanent its name will change to reflect its newly published status via the
## publication's SPIRES ID. The previous form of the analysis, possibly with
## different reference histograms, will be retained for one major version of
## Rivet with a status of "obsolete" before being removed, to give users time to
## migrate their run scripts, i.e. if an analysis is marked as obsolete in
## version 1.4.2, it will remain in Rivet's distribution until version
## 1.5.0. Obsolete analyses will not be available if Rivet is built with the
## --disable-obsolete configure flag.

lib_LTLIBRARIES += RivetALICEAnalyses.la
RivetALICEAnalyses_la_SOURCES = \
    ALICE_2010_S8624100.cc \
    ALICE_2010_S8625980.cc \
    ALICE_2010_S8706239.cc \
    ALICE_2011_S8909580.cc \
    ALICE_2011_S8945144.cc \
    ALICE_2012_I1181770.cc


lib_LTLIBRARIES += RivetATLASAnalyses.la
RivetATLASAnalyses_la_SOURCES = \
    ATLAS_2010_S8591806.cc \
    ATLAS_2010_S8817804.cc \
    ATLAS_2010_S8894728.cc \
    ATLAS_2010_S8914702.cc \
    ATLAS_2010_S8918562.cc \
    ATLAS_2010_S8919674.cc \
    ATLAS_2011_S8924791.cc \
    ATLAS_2011_S8971293.cc \
    ATLAS_2011_S8994773.cc \
    ATLAS_2011_S8983313.cc \
    ATLAS_2011_S9002537.cc \
    ATLAS_2011_S9120807.cc \
    ATLAS_2011_S9126244.cc \
    ATLAS_2011_S9128077.cc \
    ATLAS_2011_S9131140.cc \
    ATLAS_2011_S9212183.cc \
    ATLAS_2011_S9225137.cc \
    ATLAS_2011_S9019561.cc \
    ATLAS_2011_I894867.cc \
    ATLAS_2011_I919017.cc \
    ATLAS_2011_I921594.cc \
    ATLAS_2011_I925932.cc \
    ATLAS_2011_I926145.cc \
	ATLAS_2011_I930220.cc \
    ATLAS_2011_I944826.cc \
    ATLAS_2011_I945498.cc \
    ATLAS_2011_I954993.cc \
	ATLAS_2011_I928289_W.cc \
	ATLAS_2011_I928289_Z.cc \
    ATLAS_2012_I943401.cc \
    ATLAS_2012_I1082936.cc \
    ATLAS_2012_I1083318.cc \
    ATLAS_2012_I1084540.cc \
    ATLAS_2012_I1091481.cc \
    ATLAS_2012_I1093734.cc \
    ATLAS_2012_I1093738.cc \
    ATLAS_2012_I1094564.cc \
    ATLAS_2012_I1094568.cc \
    ATLAS_2012_I1112263.cc \
    ATLAS_2012_I1117704.cc \
    ATLAS_2012_I1118269.cc \
    ATLAS_2012_I1119557.cc \
    ATLAS_2012_I1124167.cc \
    ATLAS_2012_I1125575.cc \
    ATLAS_2012_I1125961.cc \
    ATLAS_2012_I1183818.cc \
    ATLAS_2012_I1188891.cc \
    ATLAS_2012_I1199269.cc \
    ATLAS_2012_I1203852.cc \
    ATLAS_2012_I1204447.cc \
    ATLAS_2012_I1204784.cc \
    ATLAS_2013_I1190187.cc \
    ATLAS_2013_I1217867.cc \
    ATLAS_2013_I1219109.cc \
    ATLAS_2013_I1230812.cc \
    ATLAS_2013_I1243871.cc \
    ATLAS_2013_I1263495.cc \
    ATLAS_2014_I1268975.cc \
	ATLAS_2014_I1279489.cc \
	ATLAS_2014_I1282441.cc \
	ATLAS_2014_I1298811.cc \
    ATLAS_2014_I1304688.cc \
    ATLAS_2014_I1307756.cc \
    ATLAS_2014_I1306294.cc \
    ATLAS_2014_I1315949.cc \
    ATLAS_2014_I1325553.cc \
    ATLAS_2014_I1300647.cc \
    ATLAS_2014_I1288706.cc \
    ATLAS_2014_I1307243.cc \
    ATLAS_2014_I1312627.cc \
	ATLAS_2013_I1217863_W.cc \
	ATLAS_2013_I1217863_Z.cc \
    ATLAS_2014_I1306615.cc
if ENABLE_PRELIMINARY
RivetATLASAnalyses_la_SOURCES += \
    ATLAS_2012_CONF_2012_001.cc
endif
if ENABLE_UNVALIDATED
RivetATLASAnalyses_la_SOURCES += \
    ATLAS_2010_CONF_2010_049.cc \
    ATLAS_2011_S9041966.cc \
    ATLAS_2011_CONF_2011_090.cc \
    ATLAS_2011_CONF_2011_098.cc \
    ATLAS_2011_S9108483.cc \
    ATLAS_2011_S9212353.cc \
    ATLAS_2011_S9035664.cc \
    ATLAS_2012_I1095236.cc \
    ATLAS_2012_I1082009.cc \
    ATLAS_2012_I946427.cc \
    ATLAS_2012_I1126136.cc \
    ATLAS_2012_I1180197.cc \
    ATLAS_2012_I1186556.cc \
    ATLAS_2012_I1190891.cc \
    ATLAS_2012_CONF_2012_103.cc \
    ATLAS_2012_CONF_2012_104.cc \
    ATLAS_2012_CONF_2012_105.cc \
    ATLAS_2012_CONF_2012_109.cc \
    ATLAS_2012_CONF_2012_153.cc
endif


lib_LTLIBRARIES += RivetCMSAnalyses.la
RivetCMSAnalyses_la_SOURCES = \
    CMS_2010_S8547297.cc \
    CMS_2010_S8656010.cc \
    CMS_2011_S8884919.cc \
    CMS_2011_S8941262.cc \
    CMS_2011_S8950903.cc \
    CMS_2011_S8957746.cc \
    CMS_2011_S8968497.cc \
    CMS_2011_S8973270.cc \
    CMS_2011_S8978280.cc \
    CMS_2011_S9086218.cc \
    CMS_2011_S9088458.cc \
    CMS_2011_S9120041.cc \
    CMS_2011_S9215166.cc \
    CMS_2012_I941555.cc  \
    CMS_2011_I954992.cc  \
    CMS_2012_I1087342.cc \
    CMS_2012_I1090423.cc \
    CMS_2012_I1102908.cc \
    CMS_2012_I1107658.cc \
    CMS_2012_I1184941.cc \
    CMS_2012_I1193338.cc \
    CMS_2013_I1209721.cc \
    CMS_2013_I1218372.cc \
    CMS_2013_I1224539_DIJET.cc \
	CMS_2013_I1224539_WJET.cc \
	CMS_2013_I1224539_ZJET.cc \
    CMS_2013_I1256943.cc \
    CMS_2013_I1258128.cc \
    CMS_2013_I1261026.cc \
    CMS_2013_I1265659.cc \
    CMS_2013_I1272853.cc \
    CMS_2013_I1273574.cc \
    CMS_2014_I1298810.cc \
    CMS_2014_I1303894.cc \
    CMSTOTEM_2014_I1294140.cc
if ENABLE_PRELIMINARY
RivetCMSAnalyses_la_SOURCES += \
    CMS_QCD_10_024.cc \
    CMS_2012_PAS_QCD_11_010.cc
endif


lib_LTLIBRARIES += RivetLHCbAnalyses.la
RivetLHCbAnalyses_la_SOURCES = \
    LHCB_2010_I867355.cc \
    LHCB_2012_I1208102.cc \
    LHCB_2013_I1218996.cc \
    LHCB_2013_I1208105.cc \
    LHCB_2014_I1281685.cc
if ENABLE_UNVALIDATED
RivetLHCbAnalyses_la_SOURCES += \
    LHCB_2010_S8758301.cc \
    LHCB_2011_I917009.cc \
    LHCB_2011_I919315.cc \
    LHCB_2012_I1119400.cc
endif


lib_LTLIBRARIES += RivetLHCfAnalyses.la
RivetLHCfAnalyses_la_SOURCES = \
    LHCF_2012_I1115479.cc


lib_LTLIBRARIES += RivetTOTEMAnalyses.la
RivetTOTEMAnalyses_la_SOURCES = \
    TOTEM_2012_I1115294.cc \
    TOTEM_2012_002.cc # TODO: update to Inspire ID


lib_LTLIBRARIES += RivetCDFAnalyses.la
RivetCDFAnalyses_la_SOURCES = \
    CDF_1988_S1865951.cc \
    CDF_1990_S2089246.cc \
    CDF_1994_S2952106.cc \
    CDF_1996_S3418421.cc \
    CDF_1998_S3618439.cc \
    CDF_2000_S4155203.cc \
    CDF_2000_S4266730.cc \
    CDF_2001_S4517016.cc \
    CDF_2001_S4563131.cc \
    CDF_2001_S4751469.cc \
    CDF_2002_S4796047.cc \
    CDF_2004_S5839831.cc \
    CDF_2005_S6080774.cc \
    CDF_2005_S6217184.cc \
    CDF_2006_S6450792.cc \
    CDF_2006_S6653332.cc \
    CDF_2007_S7057202.cc \
    CDF_2008_S7540469.cc \
    CDF_2008_S7828950.cc \
    CDF_2008_S8093652.cc \
    CDF_2008_S8095620.cc \
    CDF_2009_S8233977.cc \
    CDF_2009_S8383952.cc \
    CDF_2009_S8436959.cc \
    CDF_2010_S8591881_DY.cc \
    CDF_2010_S8591881_QCD.cc
if ENABLE_PRELIMINARY
RivetCDFAnalyses_la_SOURCES += \
    CDF_2009_NOTE_9936.cc \
    CDF_2012_NOTE10874.cc
endif
# if ENABLE_OBSOLETE
# RivetCDFAnalyses_la_SOURCES +=
# endif
if ENABLE_UNVALIDATED
RivetCDFAnalyses_la_SOURCES += \
    CDF_1993_S2742446.cc \
    CDF_1996_S3108457.cc \
    CDF_1996_S3349578.cc \
    CDF_1997_S3541940.cc \
    CDF_2008_S7541902.cc \
    CDF_2008_S7782535.cc
endif


lib_LTLIBRARIES += RivetD0Analyses.la
RivetD0Analyses_la_SOURCES = \
    D0_2000_S4480767.cc \
    D0_2000_I499943.cc \
    D0_2001_S4674421.cc \
    D0_2004_S5992206.cc \
    D0_2006_S6438750.cc \
    D0_2007_S7075677.cc \
    D0_2008_S6879055.cc \
    D0_2008_S7554427.cc \
    D0_2008_S7662670.cc \
    D0_2008_S7719523.cc \
    D0_2008_S7837160.cc \
    D0_2008_S7863608.cc \
    D0_2009_S8202443.cc \
    D0_2009_S8320160.cc \
    D0_2009_S8349509.cc \
    D0_2010_S8566488.cc \
    D0_2010_S8570965.cc \
    D0_2010_S8671338.cc \
    D0_2010_S8821313.cc \
    D0_2011_I895662.cc
if ENABLE_UNVALIDATED
RivetD0Analyses_la_SOURCES += \
    D0_1996_S3214044.cc \
    D0_1996_S3324664.cc
endif


lib_LTLIBRARIES += RivetHERAAnalyses.la
RivetHERAAnalyses_la_SOURCES = \
    H1_1994_S2919893.cc \
    H1_2000_S4129130.cc
if ENABLE_UNVALIDATED
RivetHERAAnalyses_la_SOURCES += \
    H1_1995_S3167097.cc \
    ZEUS_2001_S4815815.cc
endif

lib_LTLIBRARIES += RivetPetraAnalyses.la
RivetPetraAnalyses_la_SOURCES = \
    JADE_1998_S3612880.cc \
    TASSO_1990_S2148048.cc

lib_LTLIBRARIES += RivetLEPAnalyses.la
RivetLEPAnalyses_la_SOURCES = \
    ALEPH_1991_S2435284.cc \
    ALEPH_1996_S3486095.cc \
    ALEPH_1996_S3196992.cc \
    ALEPH_1999_S4193598.cc \
    ALEPH_2001_S4656318.cc \
    ALEPH_2002_S4823664.cc \
    ALEPH_2004_S5765862.cc \
    DELPHI_1995_S3137023.cc \
    DELPHI_1996_S3430090.cc \
    DELPHI_1999_S3960137.cc \
    DELPHI_2000_S4328825.cc \
    OPAL_1994_S2927284.cc \
    OPAL_1995_S3198391.cc \
    OPAL_1996_S3257789.cc \
    OPAL_1997_S3396100.cc \
    OPAL_1997_S3608263.cc \
    OPAL_1998_S3702294.cc \
    OPAL_1998_S3749908.cc \
    OPAL_1998_S3780481.cc \
    OPAL_2000_S4418603.cc \
    OPAL_2001_S4553896.cc \
    OPAL_2002_S5361494.cc \
    OPAL_2004_S6132243.cc \
    SLD_1996_S3398250.cc \
    SLD_1999_S3743934.cc \
    SLD_2002_S4869273.cc \
    SLD_2004_S5693039.cc
if ENABLE_PRELIMINARY
RivetLEPAnalyses_la_SOURCES += \
    DELPHI_2002_069_CONF_603.cc
endif
if ENABLE_UNVALIDATED
RivetLEPAnalyses_la_SOURCES += \
    OPAL_1993_S2692198.cc \
    DELPHI_2003_WUD_03_11.cc
endif


lib_LTLIBRARIES += RivetRHICAnalyses.la
RivetRHICAnalyses_la_SOURCES = \
    STAR_2006_S6500200.cc \
    STAR_2006_S6860818.cc \
    STAR_2006_S6870392.cc
if ENABLE_PRELIMINARY
RivetRHICAnalyses_la_SOURCES += \
    STAR_2009_UE_HELEN.cc
endif
if ENABLE_UNVALIDATED
RivetRHICAnalyses_la_SOURCES += \
    STAR_2008_S7869363.cc \
    STAR_2008_S7993412.cc
endif


lib_LTLIBRARIES += RivetSPSAnalyses.la
RivetSPSAnalyses_la_SOURCES = \
    UA1_1990_S2044935.cc \
    UA5_1982_S875503.cc \
    UA5_1986_S1583476.cc \
    UA5_1987_S1640666.cc \
    UA5_1988_S1867512.cc \
    UA5_1989_S1926373.cc


lib_LTLIBRARIES += RivetMiscAnalyses.la
RivetMiscAnalyses_la_SOURCES = \
    PDG_HADRON_MULTIPLICITIES.cc \
    PDG_HADRON_MULTIPLICITIES_RATIOS.cc \
    PDG_TAUS.cc \
    JADE_OPAL_2000_S4300807.cc \
    ARGUS_1993_S2653028.cc \
    ARGUS_1993_S2669951.cc \
    ARGUS_1993_S2789213.cc \
    BABAR_2003_I593379.cc \
    BABAR_2005_S6181155.cc \
    BABAR_2007_S6895344.cc \
    BABAR_2007_S7266081.cc \
    BABAR_2013_I1238276.cc \
    BELLE_2001_S4598261.cc \
    BELLE_2008_I786560.cc \
    BELLE_2013_I1216515.cc \
    CLEO_2004_S5809304.cc
if ENABLE_UNVALIDATED
RivetMiscAnalyses_la_SOURCES += \
    E735_1998_S3905616.cc \
    SFM_1984_S1178091.cc
endif


lib_LTLIBRARIES += RivetMCAnalyses.la
RivetMCAnalyses_la_SOURCES = \
    EXAMPLE.cc \
    EXAMPLE_CUTS.cc \
    MC_QCD_PARTONS.cc \
    MC_DIPHOTON.cc \
    MC_ELECTRONS.cc \
    MC_GENERIC.cc \
    MC_HINC.cc \
    MC_HJETS.cc \
    MC_HKTSPLITTINGS.cc \
    MC_IDENTIFIED.cc \
    MC_JETS.cc \
    MC_JETTAGS.cc \
    MC_HFJETS.cc \
    MC_LEADJETUE.cc \
    MC_MUONS.cc \
    MC_PDFS.cc \
    MC_PHOTONINC.cc \
    MC_PHOTONJETS.cc \
    MC_PHOTONKTSPLITTINGS.cc \
    MC_PHOTONS.cc \
    MC_PRINTEVENT.cc \
    MC_SUSY.cc \
    MC_TTBAR.cc \
    MC_VH2BB.cc \
    MC_WINC.cc \
    MC_WJETS.cc \
    MC_WKTSPLITTINGS.cc \
    MC_WPOL.cc \
    MC_WWINC.cc \
    MC_WWJETS.cc \
    MC_WWKTSPLITTINGS.cc \
    MC_XS.cc \
    MC_ZINC.cc \
    MC_ZJETS.cc \
    MC_ZKTSPLITTINGS.cc \
    MC_ZZINC.cc \
    MC_ZZJETS.cc \
    MC_ZZKTSPLITTINGS.cc


if ENABLE_UNVALIDATED
RivetMCAnalyses_la_SOURCES += \
    MC_DIJET.cc \
    MC_PHOTONJETUE.cc
endif
