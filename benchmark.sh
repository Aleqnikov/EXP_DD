#!/usr/bin/env bash

BINARY="./build/dd_calc"
X_VAL="${1:-124.24312412434234}"
DATA_FILE="/tmp/bench_data.dat"

echo "# iters dd_ns std_ns" > "$DATA_FILE"

for exp in $(seq 0 0.1 5); do
    n=$(python3 -c "print(max(1, int(10**$exp)))")

    out=$("$BINARY" -x "$X_VAL" -b -n "$n" 2>/dev/null)

    dd_ns=$(echo  "$out" | grep "dd::exp"  | grep -oP '[\d.e+\-]+(?= нс)')
    std_ns=$(echo "$out" | grep "std::exp" | grep -oP '[\d.e+\-]+(?= нс)')

    if [[ -n "$dd_ns" && -n "$std_ns" ]]; then
        echo "$n $dd_ns $std_ns" >> "$DATA_FILE"
    fi
done

gnuplot <<EOF
set terminal qt size 1000,600 font "Sans,12"

set title "dd::exp vs std::exp  (x = $X_VAL)"
set xlabel "Число итераций"
set ylabel "Время на вызов, нс"

set logscale x 10
set grid xtics ytics
set key top right

plot "$DATA_FILE" using 1:2 with linespoints lw 2 pt 7 ps 0.8 lc rgb "#e74c3c" title "dd::exp", \
     "$DATA_FILE" using 1:3 with linespoints lw 2 pt 7 ps 0.8 lc rgb "#2980b9" title "std::exp"

pause mouse close
EOF