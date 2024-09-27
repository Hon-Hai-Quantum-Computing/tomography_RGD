#!/bin/bash

echo " Do the measurement by each part"

Nk=3
mea=2400
StateName=GHZ
meas_path=../Tomography/DataTest/GHZ-3/GHZ-3_m36_s1/GHZ-3_m36_s1_shot2400_v2_Measure
Rec=$meas_path"/Rec_qiskit_meas_dict_"

echo "   Nk       = " $Nk
echo "   mea      = " $mea
echo " StateName  = " $StateName
echo " meas_path  = " $meas_path
echo " Rec Header = " $Rec

taskset -c 0 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_0 > $Rec"0.txt" &
taskset -c 1 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_1 > $Rec"1.txt" &
taskset -c 2 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_10 > $Rec"10.txt" &
taskset -c 3 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_11 > $Rec"11.txt" &
taskset -c 4 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_2 > $Rec"2.txt" &
taskset -c 5 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_3 > $Rec"3.txt" &
taskset -c 6 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_4 > $Rec"4.txt" &
taskset -c 7 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_5 > $Rec"5.txt" &
taskset -c 8 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_6 > $Rec"6.txt" &
taskset -c 9 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_7 > $Rec"7.txt" &
taskset -c 10 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_8 > $Rec"8.txt" &
taskset -c 11 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_9 > $Rec"9.txt" &

wait

echo "All the measurement_dict DONE"