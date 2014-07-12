#!/bin/bash
if [ ! -f ../../tpar ]; then
        ../../make
fi
if [ ! -d opt ]; then
        mkdir opt
fi
for f in `find ./ -type f | sed 's/^\.\.\///'`
do
	if [ "$f" != "GF2^128" ] && [ "$f" != "GF2^256" ] ; then
	  ../../tpar < $f > opt/$f.opt
	fi
done
