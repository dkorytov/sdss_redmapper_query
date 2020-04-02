#!/bin/bash

mstars=( -1 -0.5 0 0.5 1 )

for mstar in "${mstars[@]}"
do
    cp template.param baxter/${mstar}_crit.param
    sed -i "s/#mstar_val#/${mstar}/g" baxter/${mstar}_crit.param
    sed -i "s/#richness_mass_author#/Baxter_crit/g" baxter/${mstar}_crit.param
    sed -i "s/simet_crit/baxter_crit/g" mcclintock/${mstar}_crit.param
    sed -i "s/high_richness/mcclintock/g" mcclintock/${mstar}_crit.param
    echo $mstar
done

