grep funcT /home/magh/COIN/coin-Decomp/build-g-d/Decomp/examples/MMKP/I03.txt | awk '{if($9 > 0.01) print $0}' | sort -nk 9
