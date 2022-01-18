for i in $(seq 21 40); do
    mkdir $i
    cd $i
    ../general.x -N 256 -D 6 -s $i -n 1000000 -o a.edges > a.txt &
    cd ../
done
