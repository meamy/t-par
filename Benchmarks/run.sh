#!/bin/bash
if [ ! -f ../t-par ]; then
        ../make
fi
if [ ! -d opt ]; then
        mkdir opt
fi
for f in `find ./*.qc`
do
  ../t-par < $f > opt/$f.opt
done
