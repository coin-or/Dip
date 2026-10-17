cd ../../ThirdParty

wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/hmetis/hmetis-1.5-linux.tar.gz

gunzip hmetis-1.5-linux.tar.gz

tar -xvf hmetis-1.5-linux.tar

rm hmetis-1.5-linux.tar

mv hmetis-1.5-linux HMetis

cd HMetis

wget http://coral.ie.lehigh.edu/~jiadongwang/DIP/hmetis.h
