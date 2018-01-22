#!/bin/sh
fiso=A549
fcon=xglc
inpar="glc/1"
oudir="glc/"
fstat="glc/statfl"
fcmpr="glut/statfl"
manfi=77
fit=false
numpar=false
confint=false
tkfit=false
tafit=false
dostat=false
tst=yes

while getopts ":a:b:i:o:s:c:m:F:N:C:K:A:S:" opt; do
  case $opt in
    a) fiso=$OPTARG;;
    b) fcon=$OPTARG;;
    i) inpar=$OPTARG;;
    o) oudir=$OPTARG;;
    s) fstat=$OPTARG;;
    c) fcmpr=$OPTARG;;
    m) manfi=$OPTARG;;
    F) fit=$OPTARG;;
    N) numpar=$OPTARG;;
    C) confint=$OPTARG;;
    K) tkfit=$OPTARG;;
    A) tafit=$OPTARG;;
    S) dostat=$OPTARG;;
    *)
      echo "Invalid option: -$OPTARG" 
      cat help
      tst=no
      ;;
  esac
done
if [ $tst = yes ]
then
./isodyn.out $fiso $fcon $inpar $oudir $fstat $fcmpr $manfi $fit $numpar $confint $tkfit $tafit $dostat
fi
