#!/usr/bin/env bash

VER="3.1.2"

# use Python version from Docker image
alias python='docker run -i --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:'$VER' python'

# Rivet tools ...
alias rivet='docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:'$VER' rivet'
alias rivet-build='docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:'$VER' rivet-build'
alias rivet-cmphistos='docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:'$VER' rivet-cmphistos'
alias rivet-mkanalysis='docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:'$VER' rivet-mkanalysis'
alias rivet-mkhtml='docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:'$VER' rivet-mkhtml'
alias make-plots='docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:'$VER' make-plots'

# YODA tools ...
alias yodamerge='docker run -i --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:'$VER' yodamerge'
alias yodascale='docker run -i --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:'$VER' yodascale'
alias yodadiff='docker run -i --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:'$VER' yodadiff'
alias yoda2flat='docker run -i --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:'$VER' yoda2flat'
alias yodals='docker run -i --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:'$VER' yodals'

# let Rivet find local libraties
export RIVET_ANALYSIS_PATH=$PWD

