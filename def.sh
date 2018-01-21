#!/bin/bash

while getopts ":a" opt; do
  case $opt in
    a)
      echo "-a was triggered!" 
shift $((OPTIND-1))
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done
variable1=$1
variable2=${2:-asd}
echo $variable1
echo $variable2

