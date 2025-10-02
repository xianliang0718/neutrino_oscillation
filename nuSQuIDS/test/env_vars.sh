
export CXX="g++"
export CFLAGS=" -I/home/xianliang/download/musairs/nuSQuIDS/inlude -I/usr/local/include -Wno-abi   -I/usr/include/hdf5/serial"
export CXXFLAGS=" -std=c++11"
export LDFLAGS=" -L/home/xianliang/download/musairs/nuSQuIDS/lib -Wl,-rpath -Wl,/home/xianliang/download/musairs/nuSQuIDS/lib -lnuSQuIDS -lpthread -L/usr/local/lib -lSQuIDS -lgsl -lgslcblas -lm  -lgsl -lgslcblas -lm  -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_hl -lhdf5 -lcrypto -lcurl -lpthread -lsz -lz -ldl -lm"

export LD_LIBRARY_PATH="/home/xianliang/download/musairs/nuSQuIDS/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu/hdf5/serial:/home/xianliang/download/globes/globes-3.2.18//lib:/opt/geant4/lib:/opt/root/lib:/home/xianliang/download/globes/globes-3.2.18//lib:/opt/geant4/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu"
