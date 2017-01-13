MPICC=mpicc

default:
	$(MPICC) -std=c99 -O3 -march=native src/main.c -o many_body -lm