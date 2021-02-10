DRAW_PY=../../scripts/draw_graph.py
TIMES=5
rm -f *.edges *.png
gcc grid.c ../../src/libodp.a -I.. -o grid.x
echo 1
for i in $(seq 1 $TIMES); do
    ./grid.x -w 6 -h 6 -d 3 -l 3 -g 1 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
for i in $(seq 1 $TIMES); do
    ./grid.x -w 10 -h 10 -d 4 -l 4 -g 1 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
echo 2
rm -f *.edges *.png
for i in $(seq 1 $TIMES); do
    ./grid.x -w 6 -h 6 -d 3 -l 3 -g 2 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
for i in $(seq 1 $TIMES); do
    ./grid.x -w 8 -h 8 -d 4 -l 3 -g 2 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
echo 4
rm -f *.edges *.png
for i in $(seq 1 $TIMES); do
    ./grid.x -w 6 -h 6 -d 3 -l 3 -g 4 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
for i in $(seq 1 $TIMES); do
    ./grid.x -w 8 -h 8 -d 4 -l 3 -g 4 -s $i > $i.edges; python $DRAW_PY $i.edges -n
done
rm -f *.edges *.png grid.x
