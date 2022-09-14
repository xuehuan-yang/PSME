#!/bin/bash
if [ $# -eq 0 ]
then 
    pushmessage=`date +%Y-%m-%d@%H:%M:%S-%A`
else
    pushmessage="$*"
fi

echo ${pushmessage}

git pull;
git status;
git add -A; 
git commit -m "${pushmessage}";
git push -u origin main;
