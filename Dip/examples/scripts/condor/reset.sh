rm -f ${HOME}/running/decomp${1}
mkdir ${HOME}/running/decomp${1}
for i in atm gap milpblock mmkp
do
  mkdir ${HOME}/running/decomp${1}/${i}
  mkdir ${HOME}/running/decomp${1}/${i}/${i}p
  mkdir ${HOME}/running/decomp${1}/${i}/${i}p2
  mkdir ${HOME}/running/decomp${1}/${i}/${i}c
  mkdir ${HOME}/running/decomp${1}/${i}/${i}c2
done
