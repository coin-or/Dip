cd ../../ThirdParty

wget http://bmi.osu.edu/~umit/PaToH/patoh-Linux-x86_64.tar.gz

gunzip patoh-Linux-x86_64.tar.gz

mkdir PaToH_v3.2

tar xvf patoh-Linux-x86_64.tar -C PaToH_v3.2/

rm *.tar

cd ./PaToH_v3.2/

cp -rf ./build/* .

rm -rf ./build
