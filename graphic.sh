#!/usr/bin/env bash

BINARY="./build/dd_calc"
START="${1:-0}"
END="${2:-10}"
STEP="${3:-0.1}"
DATA_FILE="/tmp/dd_exp_plot.dat"

ADD_STD=0
if [[ "$*" =~ -c ]]; then
    ADD_STD=1
fi

if [[ "$ADD_STD" -eq 1 ]]; then
    echo "# x dd std" > "$DATA_FILE"
else
    echo "# x dd" > "$DATA_FILE"
fi

awk -v start="$START" -v end="$END" -v step="$STEP" '
BEGIN {
    for (x = start; x <= end; x += step)
        printf "%.10g\n", x
}
' | while read -r x; do

    out=$("$BINARY" -x "$x" $([[ "$ADD_STD" -eq 1 ]] && echo "-c") 2>/dev/null)

    hi=$(echo "$out" | awk '/hi/ {for(i=1;i<=NF;i++) if($i=="=") print $(i+1)}')
    lo=$(echo "$out" | awk '/lo/ {for(i=1;i<=NF;i++) if($i=="=") print $(i+1)}')

    if [[ -n "$hi" && -n "$lo" ]]; then
        val=$(awk "BEGIN {print $hi + $lo}")

        if [[ "$ADD_STD" -eq 1 ]]; then
            std=$(echo "$out" | awk '/std/ {print $NF}')
            if [[ -n "$std" ]]; then
                echo "$x $val $std" >> "$DATA_FILE"
            fi
        else
            echo "$x $val" >> "$DATA_FILE"
        fi
    fi

done


if [[ "$ADD_STD" -eq 1 ]]; then
gnuplot <<EOF
set terminal qt size 1000,600
set title "exp(x)"
plot "$DATA_FILE" using 1:2 with lines title "dd", \
     "$DATA_FILE" using 1:3 with lines title "std"
pause mouse close
EOF
else
gnuplot <<EOF
set terminal qt size 1000,600
set title "exp(x)"
plot "$DATA_FILE" using 1:2 with lines title "dd"
pause mouse close
EOF
fi