compile:
	gcc -o addSolvent addSolvent.c -lm -Wall -O3
all:
	gcc -o addSolvent addSolvent.c -lm -Wall -O3
	./addSolvent dump.lammpstrj output.data solvated
help:
	./addSolvent
