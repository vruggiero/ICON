#!/bin/ksh

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

# -------------------------------------------------
# ICON atmospheric water  budget
#
# run as: water_budget.s
#
# dM/dt = E + P [kg/m2/s]   M: atm. water mass, E: surface evaporation, P: surface precip
#
# Info: create output files
#       - file_flux: fluxes time-mean between output steps with operation=mean
#                    rain_gsp_rate, snow_gsp_rate, ice_gsp_rate, (graupel_gsp_rate),
#                    rain_con_rate, snow_con_rate, qhfl_s, qcfl_s, qifl_s, snow_blow
#       - file_stat: state variables: tqv,tqc,tqi,tqr,tqs (tqg)
#       - output_grid=.TRUE.  (add grid file for weights in fldmean)
# History of errors (MAE):
#       - 21.4.2024 (f629e1ff8):  mae=3.1226e-12
#
# Martin Koehler, Roland Wirth, 4-2024
# -------------------------------------------------

#datadir=/hpc/uwork/mkoehler/run-icon/experiments/couple_090
#datadir=${EXPDIR}

for irun in 0 2; do

    tqg=''
    graupel_gsp_rate=''

    file_flux=water_budget_flux_${irun}_DOM01_ML_0001.nc
    file_stat=water_budget_state_${irun}_DOM01_ML_0001.nc

    dtime=360  #86400      # output interval [s]

    if [[ irun == 0 ]]; then
        threshold="1*10^-11"     # threshold in MAE of water budget
    else
        threshold="2*10^-11"
    fi

    # -------------------------------------------------

    cdo="${cdo:-/hpc/sw/cdo/2.4.0/x86/gnu/bin/cdo}"

    echo "=== CHECK RUN ${irun} ==="
    echo '--- input files: ' ${file_flux} ${file_stat}

    #cd ${datadir}

    #-- difference between time steps for all water integrals (tqv, tqc, tqi, tqs)
    $cdo deltat ${file_stat} tqx.nc

    #-- total atmospheric water (tqtot) and its change (dt_tqtot)
    $cdo -fldmean -expr,"dt_tqtot = (tqv + tqc + tqi + tqr + tqs $tqg) / ${dtime}"  tqx.nc    dt_tqtot.nc

    #-- Precip-Evap including turbulent qc flux at surface (qcfl_s)
    $cdo -fldmean -expr,"pme = rain_con_rate_corr + snow_con_rate_corr + ice_gsp_rate + snow_gsp_rate + rain_gsp_rate $graupel_gsp_rate + qhfl_s + qcfl_s" \
         -delete,timestep=1 ${file_flux} pme.nc

    #-- total water budget: ratio and total
    $cdo -O merge dt_tqtot.nc pme.nc temp.nc
    $cdo -expr,'ratio = - pme / dt_tqtot' temp.nc temp2.nc
    $cdo -expr,'total =   pme + dt_tqtot' temp.nc temp3.nc

    $cdo -O merge dt_tqtot.nc pme.nc temp2.nc temp3.nc budget.nc

    #-- mean absolute error (MSE) in water budget (time mean of absolute error - total)
    $cdo timmean -abs budget.nc budget_mean.nc

    echo ''
    echo '--- atmospheric water budget [kg/m2/s] --- :'
    echo '    pme:      precip-evap'
    echo '    dt_tqtot: atm. storge change in output step'
    echo '    ratio:  - pme / dt_tqtot'
    echo '    total:    pme + dt_tqtot'

    #$cdo infon pme.nc
    #$cdo infon dtqtot.nc
    #$cdo infon all2.nc
    $cdo infon budget.nc | grep total
    $cdo infon budget_mean.nc

    # --- clean-up
    rm -f tqx.nc temp.nc temp2.nc temp3.nc # dt_tqtot.nc pme.nc


    # --- test if mean absolut error (total) is smaller than a threshold
    # --- total = pme + dt_in kg/m2/s is smaller

    mae=$($cdo infon budget_mean.nc | grep 'total' | awk '{print $9}')

    echo MAE: ${mae}   threshold: ${threshold}

    # Replace 'e' with '*10^' in the value and threshold
    mae_2=$(echo $mae | sed 's/[eE]/*10^/')

    result=$(echo "$mae_2 > $threshold" | bc -l)

    if (( result != 0 )); then
        echo "MAE of run ${irun} is larger than $threshold."
        exit 1
    fi
done

echo "Success: MAE is less than $threshold."
