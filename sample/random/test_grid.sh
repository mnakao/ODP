DRAW_PY=../../script/draw_graph.py
TIMES=5
rm -f *.edges *.png
gcc grid.c ../../libodp.a -I../../include -o grid.x
echo 1
for i in $(seq 1 $TIMES); do
    ./grid.x -w 6 -h 6 -d 3 -l 3 -g 1 -s $i > $i.edges; python3 $DRAW_PY $i.edges -n
done
sleep 1
for i in $(seq 1 $TIMES); do
    ./grid.x -w 10 -h 10 -d 4 -l 4 -g 1 -s $i > $i.edges; python3 $DRAW_PY $i.edges -n
done
sleep 1
echo 2
rm -f *.edges *.png
for i in $(seq 1 $TIMES); do
    ./grid.x -w 6 -h 6 -d 3 -l 3 -g 2 -s $i > $i.edges; python3 $DRAW_PY $i.edges -n
done
sleep 1
for i in $(seq 1 $TIMES); do
    ./grid.x -w 8 -h 8 -d 4 -l 3 -g 2 -s $i > $i.edges; python3 $DRAW_PY $i.edges -n
done
sleep 1
echo 4
rm -f *.edges *.png
for i in $(seq 1 $TIMES); do
    ./grid.x -w 6 -h 6 -d 3 -l 3 -g 4 -s $i > $i.edges; python3 $DRAW_PY $i.edges -n
done
sleep 1
for i in $(seq 1 $TIMES); do
    ./grid.x -w 8 -h 8 -d 4 -l 3 -g 4 -s $i > $i.edges; python3 $DRAW_PY $i.edges -n
done
rm -f *.edges *.png grid.x
