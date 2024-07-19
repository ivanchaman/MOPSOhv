CC=gcc
CFLAGS=
LDFLAGS=
OBJS= 

all:            mopsocd

mopsocd:	mopsocd.c test-fun.h
		$(CC) $(CFLAGS) $(LDFLAGS) -o mopsocd mopsocd.c $(OBJS) -lm

clean:
		rm *.o 

real-clean:          
		rm *~ *.o *.out
