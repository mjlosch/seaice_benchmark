#!/bin/bash

dstr='--d T,172800.0'
fname='snapshot.0000000000.t001.nc'

rname='runfe01'
extns='1440'
ncks -d T,172800.0 2000/${rname}/${fname} mitgcm2km_${extns}.cdf
ncks -d T,172800.0 4000/${rname}/${fname} mitgcm4km_${extns}.cdf
ncks -d T,172800.0 8000/${rname}/${fname} mitgcm8km_${extns}.cdf
ncks -d T,172800.0 16000/${rname}/${fname} mitgcm16km_${extns}.cdf

rname='runfe02'
extns='1440_adv7'
ncks -d T,172800.0 2000/${rname}/${fname} mitgcm2km_${extns}.cdf
ncks -d T,172800.0 4000/${rname}/${fname} mitgcm4km_${extns}.cdf
ncks -d T,172800.0 8000/${rname}/${fname} mitgcm8km_${extns}.cdf
ncks -d T,172800.0 16000/${rname}/${fname} mitgcm16km_${extns}.cdf

# ncrcat mitgcm16km_1440*.cdf 16.cdf
# ncrcat mitgcm8km_1440*.cdf 8.cdf
# ncrcat mitgcm4km_1440*.cdf 4.cdf
# ncrcat mitgcm2km_1440*.cdf 2.cdf

rname='runfe03'
extns='1440_adv1'
ncks -d T,172800.0 2000/${rname}/${fname} mitgcm2km_${extns}.cdf

tar cf mitgcm_output.tar mitgcm*.cdf
