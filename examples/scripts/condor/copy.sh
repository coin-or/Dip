OPT=$1       #example=-g or -O
VERSION=$2   #example=   or -10 (the latter uses CPA10.2)

echo "Copy MMKP ${OPT}"
cp ../../../../build${OPT}${VERSION}/Dip/examples/MMKP/decomp_mmkp ${HOME}/bin${VERSION}/decomp/decomp_mmkp${OPT}
echo "Copy ATM ${OPT}"
cp ../../../../build${OPT}${VERSION}/Dip/examples/ATM/decomp_atm ${HOME}/bin${VERSION}/decomp/decomp_atm${OPT}
echo "Copy GAP ${OPT}"
cp ../../../../build${OPT}${VERSION}/Dip/examples/GAP/decomp_gap ${HOME}/bin${VERSION}/decomp/decomp_gap${OPT}
echo "Copy MilpBlock ${OPT}"
cp ../../../../build${OPT}${VERSION}/Dip/examples/MILPBlock/decomp_milpblock ${HOME}/bin${VERSION}/decomp/decomp_milpblock${OPT}
