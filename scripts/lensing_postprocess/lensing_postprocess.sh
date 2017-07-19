#!/bin/bash

START=7
END=7

for i in $(seq $START $END);
do
  DIR="k17_S${i}"
  
  mkdir -p $DIR
  
  mv "dust_lensing_R-128_F-800_kcut-17_Nside-16_S-${i}" $DIR
  mv "dust_lensing_R-160_F-1000_kcut-17_Nside-16_S-${i}" $DIR
  mv "dust_lensing_R-192_F-1200_kcut-17_Nside-16_S-${i}" $DIR
  mv "dust_lensing_R-256_F-1600_kcut-17_Nside-16_S-${i}" $DIR

  cp raytrace2k.m $DIR

  cd $DIR
  math < raytrace2k.m
  cd ..

  cp "k17_S${i}/grf_Nside16_kappx.txt" "txtfiles/grf_Nside16_k17_S${i}_kappx.txt"
  cp "k17_S${i}/grf_Nside16_kGR.txt" "txtfiles/grf_Nside16_k17_S${i}_kGR.txt"
  cp "k17_S${i}/grf_Nside16_kGRminuskappx.txt" "txtfiles/grf_Nside16_k17_S${i}_kGRminuskappx.txt"

done

