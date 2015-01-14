NAME: Timescale Distribution and Mass Measurements of Microlensing Pulsar Events
AUTHOR: SHI DAI
VERSION: 2.0 2-MAR-2014

UPDATES: New scale height of PSR distribution has been used (Lorimer et al. 2006)

mpicc -o timescale_all.out timescale_all.c -Wall -I/usr/local/include -L/usr/local/lib -lm -L/psr/dist/x86_64_lenny/lib/ -I/psr/dist/x86_64_lenny/include/ -lgsl -lgslcblas

mpirun -n 4 ./timescale_all.out
