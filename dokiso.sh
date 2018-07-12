#!/bin/sh
fiso=/data/SW620-Glucose 
fcon=/data/xglc
inpar="/data/glc/1"
oudir="/data/glc"
fstat="/data/glc/statfl"
fcmpr="/data/glut/statfl"
manfi=2
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
sudo docker run -i -t -v $PWD:/data isodyn:0.2 $fiso $fcon $inpar $oudir $fstat $fcmpr $manfi $FNCKAS
fi
