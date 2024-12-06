tar -xjf htslib-1.21.tar.bz2
cd htslib-1.21
./configure --prefix=`pwd`
make
make install
cd ..

unzip sparsehash-sparsehash-2.0.3.zip
cd sparsehash-sparsehash-2.0.3
./configure --prefix=`pwd`
make
make install
