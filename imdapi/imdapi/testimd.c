/*
 * A simple example program that implements the IMD protocol 
 * used by VMD and NAMD.
 *
 */
#include "vmdsock.h"
#include "imd.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main (int argc, char **argv) {
  int port = 54321;
  void *sock;
  void *clientsock;
  int length;
  IMDEnergies energies;
  float coords[3];
  float tmp;

  printf("%s: setting up incoming socket\n", argv[0]);
  vmdsock_init();
  sock = vmdsock_create();
  clientsock = NULL;
  vmdsock_bind(sock, port);

  printf("%s: waiting for IMD connection on port %d...\n", argv[0], port);
  vmdsock_listen(sock);

  while (!clientsock) {
    if (vmdsock_selread(sock, 0) > 0) {
      clientsock = vmdsock_accept(sock);
      if (imd_handshake(clientsock)) {
        clientsock = NULL;
      };
    }
  }

  sleep(1);
  if (vmdsock_selread(clientsock, 0) != 1 ||
      imd_recv_header(clientsock, &length) != IMD_GO) {
    clientsock = NULL;
  }

  /* jitter atom coordinate until the connection is terminated */
  tmp = 0.23234; 
  while (clientsock) {
    imd_send_energies(clientsock, &energies);
    tmp = ( (double)rand() / (double)(RAND_MAX) );
    coords[0] = tmp ;			/* coordinate x */
    coords[1] = (tmp -1);		/* coordinate y */
    coords[2] = (-tmp + 0.8);		/* coordinate z */
    imd_send_fcoords(clientsock, 1, coords);
  }
}

