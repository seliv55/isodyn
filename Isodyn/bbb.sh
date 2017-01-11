#!/bin/bash
pkill a549.out
rm nohup.out *~ */*~
./a549.out .
./limpiar.out 
