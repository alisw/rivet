
# define custom colours
DBLUE="blue!90!white"
MBLUE="blue!60!white"
LBLUE="blue!30!white"
# some useful style options for bands from https://gitlab.com/hepcedar/rivet/-/blob/release-3-1-x/doc/tutorials/makeplots.md
COMMON_TAGS="RatioPlotSameStyle=1:ErrorBars=0:ErrorBands=1:ErrorBandOpacity=0.5"
MEPS_STYLE="$COMMON_TAGS:ErrorBandColor=$LBLUE:LineColor=$LBLUE"
ME_STYLE="$COMMON_TAGS:ErrorBandColor=$MBLUE:LineColor=$MBLUE"
STAT_STYLE="$COMMON_TAGS:ErrorBandStyle=hlines:ErrorBandColor=$DBLUE:LineColor=$DBLUE"
######################################################################
# Complete regular expressions here to select the correct muliweights
######################################################################
MEPS_VARS="^MU.*PDF.*"
ME_VARS="ME_ONLY_MU.*PDF.*"
######################################################################
MEPS_BAND="Variations=$MEPS_VARS:BandComponentEnv=$MEPS_VARS"
ME_BAND="Variations=$ME_VARS:BandComponentEnv=$ME_VARS"
NO_BAND="Variations=none"
# use "Name" tag declare multiple distinct curves using the same file,
# apply style options, set a title, select the multiweights to be 
# combined into bands
rivet-mkhtml --errs -o plots_with_bands3 \
Rivet.yoda":Name=full:$MEPS_STYLE:$MEPS_BAND:Title=ME+PS scales \$\\oplus\$ stats" \
Rivet.yoda":Name=meonly:$ME_STYLE:$ME_BAND:Title=ME scales \$\\oplus\$ stats" \
Rivet.yoda":Name=stat:$STAT_STYLE:$NO_BAND:Title=stats only"

