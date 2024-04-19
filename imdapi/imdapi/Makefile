##
## Makefile for Simple Interactive MD interface example program
##

CC = cc
CXX = CC

testimd : testimd.c vmdsock.c imd.h imd.c
	$(CC) $(CFLAGS) testimd.c vmdsock.c imd.c -o testimd  -lnsl -lsocket
	@echo "See the README file for info on running this test code..."

clean :
	rm -f *.o testimd
