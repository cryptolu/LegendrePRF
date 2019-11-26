# example:
# make threads=8 target=P40
# make threads=8 target=P64
# make threads=24 target=P74

threads = 8
target = P40
CXX=clang++

buildrun:
	$(CXX) -D$(target) -DNTHREADS=$(threads) -std=c++11 -lstdc++ solve.cpp -o solve -Ofast -lpthread -lgmp
	time ./solve
