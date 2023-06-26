PROC = 8
THREADS = 8

build:
	mpic++ -DNUM_THREADS=$(THREADS) -std=c++11 -o lca_distributed lca_distributed.cpp -fopenmp

run:
	mpirun -np $(PROC) ./lca_distributed

clean:
	rm lca_distributed
