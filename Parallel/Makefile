THREADS = 8
SCHEDULE = auto

build:
	g++ -g -DSCHEDULE=$(SCHEDULE) -DNUM_THREADS=$(THREADS) lca_parallel.cpp -o lca_parallel -fopenmp

run:
	./lca_parallel

clean:
	rm lca_parallel
