line=$(head -n 1 dates.txt)

yearstr=${line:0:4}
monthstr=${line:4:2}
daystr=${line:6:2}
cyclestr=${line:8:2}

echo $yearstr
echo $monthstr
echo $daystr
echo $cyclestr

# Remove top line from dates file
tail -n +2 dates.txt > dates2.txt
mv dates2.txt dates.txt
