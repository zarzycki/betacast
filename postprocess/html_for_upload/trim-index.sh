#!/bin/bash

FIRSTLINETODELETE=77

#####

sed "${FIRSTLINETODELETE},\$d" index.HOLD > index.tmp
echo >> index.tmp
echo "</body>" >> index.tmp
echo "</html>" >> index.tmp

mv index.tmp index.HOLD



