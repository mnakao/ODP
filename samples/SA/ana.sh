for s in 1 2 3 4 6 8 9 12; do
    grep "ASPL Gap" s$s.txt | grep -v Ncalcs | awk '{print $4}' > $s
done

paste 1 2 3 4 6 8 9 12
rm -f 1 2 3 4 6 8 9 12
