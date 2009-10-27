rm -f ${HOME}/running/decomp${1}
mkdir ${HOME}/running/decomp${1}
for i in atm gap milpblock mmkp
do
  mkdir ${HOME}/running/decomp${1}/${i}p
  mkdir ${HOME}/running/decomp${1}/${i}p2
  mkdir ${HOME}/running/decomp${1}/${i}c
  mkdir ${HOME}/running/decomp${1}/${i}c2
done
