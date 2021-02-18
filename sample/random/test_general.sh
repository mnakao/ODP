DRAW_PY=../../scripts/draw_graph.py
TIMES=5
rm -f *.edges *.png
gcc general.c ../../src/libodp.a -I.. -o general.x

echo 1
for i in $(seq 1 $TIMES); do
    ./general.x -n 40 -d 3 -g 1 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
echo 2
rm -f *.edges *.png
for i in $(seq 1 $TIMES); do
    ./general.x -n 30 -d 3 -g 2 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
echo 3
rm -f *.edges *.png
for i in $(seq 1 $TIMES); do
    ./general.x -n 20 -d 3 -g 3 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
echo 4
rm -f *.edges *.png
for i in $(seq 1 $TIMES); do
    ./general.x -n 18 -d 3 -g 4 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
echo 5
rm -f *.edges *.png
for i in $(seq 1 $TIMES); do
    ./general.x -n 12 -d 4 -g 5 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
echo 6
rm -f *.edges *.png
for i in $(seq 1 $TIMES); do
    ./general.x -n 12 -d 3 -g 6 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
echo 7
rm -f *.edges *.png
for i in $(seq 1 $TIMES); do
    ./general.x -n 10 -d 3 -g 7 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
echo 8
rm -f *.edges *.png
for i in $(seq 1 $TIMES); do
    ./general.x -n 8 -d 3 -g 8 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
echo 9
rm -f *.edges *.png
for i in $(seq 1 $TIMES); do
    ./general.x -n 8 -d 3 -g 9 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
echo 10
rm -f *.edges *.png
for i in $(seq 1 $TIMES); do
    ./general.x -n 7 -d 4 -g 10 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
echo 11
rm -f *.edges *.png
for i in $(seq 1 $TIMES); do
    ./general.x -n 6 -d 4 -g 11 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
rm -f *.edges *.png general.x
