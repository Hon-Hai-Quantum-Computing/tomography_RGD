#!/bin/bash

echo " Do the measurement by each part"

Nk=4
mea=4800
StateName=GHZ
meas_path=calc/GHZ-4/GHZ-4_m256_s6/GHZ-4_m256_s6_shot4800_v1_Measure
Rec=$meas_path"/Rec_qiskit_meas_dict_"

echo "   Nk       = " $Nk
echo "   mea      = " $mea
echo " StateName  = " $StateName
echo " meas_path  = " $meas_path
echo " Rec Header = " $Rec

taskset -c 0 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_5 > $Rec"5.txt" &
taskset -c 1 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_15 > $Rec"15.txt" &
taskset -c 2 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_26 > $Rec"26.txt" &
taskset -c 3 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_14 > $Rec"14.txt" &
taskset -c 4 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_2 > $Rec"2.txt" &
taskset -c 5 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_28 > $Rec"28.txt" &
taskset -c 6 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_18 > $Rec"18.txt" &
taskset -c 7 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_7 > $Rec"7.txt" &
taskset -c 8 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_12 > $Rec"12.txt" &
taskset -c 9 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_22 > $Rec"22.txt" &
taskset -c 10 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_20 > $Rec"20.txt" &
taskset -c 11 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_4 > $Rec"4.txt" &
taskset -c 12 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_25 > $Rec"25.txt" &
taskset -c 13 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_24 > $Rec"24.txt" &
taskset -c 14 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_29 > $Rec"29.txt" &
taskset -c 15 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_1 > $Rec"1.txt" &
taskset -c 16 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_31 > $Rec"31.txt" &
taskset -c 17 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_19 > $Rec"19.txt" &
taskset -c 18 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_0 > $Rec"0.txt" &
taskset -c 19 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_13 > $Rec"13.txt" &
taskset -c 20 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_3 > $Rec"3.txt" &
taskset -c 21 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_30 > $Rec"30.txt" &
taskset -c 22 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_6 > $Rec"6.txt" &
taskset -c 23 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_10 > $Rec"10.txt" &
taskset -c 24 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_21 > $Rec"21.txt" &
taskset -c 25 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_11 > $Rec"11.txt" &
taskset -c 26 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_27 > $Rec"27.txt" &
taskset -c 27 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_9 > $Rec"9.txt" &
taskset -c 28 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_8 > $Rec"8.txt" &
taskset -c 29 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_16 > $Rec"16.txt" &
taskset -c 30 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_23 > $Rec"23.txt" &
taskset -c 31 python parallel_measurements.py $Nk $mea $StateName $meas_path ml_labels_17 > $Rec"17.txt" &

wait

echo "All the measurement_dict DONE"