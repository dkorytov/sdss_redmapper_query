#!/bin/bash

mstars=( -1 -0.5 0 0.5 1 )

for mstar in "${mstars[@]}"
do
    cp template.param high_richness/${mstar}_simet_crit.param
    sed -i "s/#mstar_val#/${mstar}/g" high_richness/${mstar}_simet_crit.param
    sed -i "s/#richness_mass_author#/Simet_crit/g" high_richness/${mstar}_simet_crit.param
    echo $mstar
done

