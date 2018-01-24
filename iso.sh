#!/bin/sh
fiso=A549
fcon=xglc
inpar="glc/1"
oudir="glc/"
fstat="glc/statfl"
fcmpr="glut/statfl"
manfi=77
FNCKAS="0"
tst=yes

while getopts ":a:b:i:o:s:c:m:FNCKAS" opt; do
  case $opt in
    a) fiso=$OPTARG;;
    b) fcon=$OPTARG;;
    i) inpar=$OPTARG;;
    o) oudir=$OPTARG;;
    s) fstat=$OPTARG;;
    c) fcmpr=$OPTARG;;
    m) manfi=$OPTARG;;
    F) FNCKAS=F;;
    N) FNCKAS=N;;
    C) FNCKAS=C;;
    K) FNCKAS=K;;
    A) FNCKAS=A;;
    S) FNCKAS=S;;
    *)
      echo "Invalid option: -$OPTARG" 
      cat help
      tst=no
      ;;
  esac
done
if [ $tst = yes ]
then
./isodyn.out $fiso $fcon $inpar $oudir $fstat $fcmpr $manfi $FNCKAS
fi
