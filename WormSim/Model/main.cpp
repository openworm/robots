// WormSim main.

#include <stdio.h>
#include "wormsim.h"

int main(int argc, char *argv[])
{
	if (!init_test())
	{
	  fprintf(stderr, "Cannot initialize\n");
	  return 1;
	}
	while (step()) {}
	term();
	return 0;
}
