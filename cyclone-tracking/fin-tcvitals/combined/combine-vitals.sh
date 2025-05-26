#!/bin/bash

cat combined_tcvitals.????.dat > ALL_combined_tcvitals.dat

gzip -f ALL_combined_tcvitals.dat
