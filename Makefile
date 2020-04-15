CC = gcc -fopenmp
LIBS = -lm
all: commands.h main.o functions.o
	$(CC) main.o functions.o -o 2loop $(LIBS)

functions.o: functions.c
	$(CC) -c functions.c $(LIBS)

main.o: main.c
	$(CC) -c main.c $(LIBS)
clean:
	rm -f *.o 2loop ./output/*.txt
