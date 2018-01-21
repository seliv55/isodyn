#!/bin/sh
fiso=A549
fcon=xglc
inpar="glc/1"
oudir="glc/"
fstat="glc/statfl"
fcmpr="glut/statfl"
manfi=77


while getopts ":a:b:i:o:s:c:m:" opt; do
  case $opt in
    a) fiso=$OPTARG;;
    b) fcon=$OPTARG;;
    i) inpar=$OPTARG;;
    o) oudir=$OPTARG;;
    s) fstat=$OPTARG;;
    c) fcmpr=$OPTARG;;
    m) manfi=$OPTARG;;
    \?)
      echo "Invalid option: -$OPTARG" 
      ;;
  esac
done
./isodyn.out $fiso $fcon $inpar $oudir $fstat $fcmpr 

