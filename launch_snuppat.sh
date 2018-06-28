# Warning
say="Press return to continue. If there is an output/processed.dat file it "
say+="will be erased if the simulation starts from the beginning."

# Warn and pause
echo $say && read -p ""

if [ $# -eq 0 ]; then
    ./snuppat
else

    mpirun -n $1 snuppat
fi
