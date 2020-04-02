#!/bin/bash

mstars=( -1 -0.5 0 0.5 1 )

for mstar in "${mstars[@]}"
do
    cp template.param farahi/${mstar}_crit.param
    sed -i "s/#mstar_val#/${mstar}/g" farahi/${mstar}_crit.param
    sed -i "s/#richness_mass_author#/Farahi_crit/g" farahi/${mstar}_crit.param
    sed -i "s/simet_crit/farahi_crit/g" farahi/${mstar}_crit.param
    sed -i "s/high_richness/farahi/g" farahi/${mstar}_crit.param
    echo $mstar
done

