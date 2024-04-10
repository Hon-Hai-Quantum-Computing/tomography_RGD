#!/bin/bash

echo " Do the measurement by each part"

Nk=3
mea=2400
StateName=GHZ
meas_path=../Tomography/DataTest/GHZ-3/GHZ-3_m10_s1/GHZ-3_m10_s1_shot2400_v1_Measure
Rec=$meas_path"/Rec_qiskit_meas_dict_"

echo "   Nk       = " $Nk
echo "   mea      = " $mea
echo " StateName  = " $StateName
echo " meas_path  = " $meas_path
echo " Rec Header = " $Rec

taskset -c 0 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_2 > $Rec"2.txt" &
taskset -c 1 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_1 > $Rec"1.txt" &
taskset -c 2 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_0 > $Rec"0.txt" &

wait

echo "All the measurement_dict DONE"