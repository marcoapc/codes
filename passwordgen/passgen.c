#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void Instructions(void);

int main(int argc, char *argv[]) {
	int i;
	int typ; //0-symbol 1-number 2-lowercase 3-uppercase
	int ntype = 4;
	int initial[4] = {33, 48, 97, 65};
	int tsize[4] = {6, 10, 26, 26}; //6
	int num_el, //number of elements in the password
	    r; //generated random number

	// Testing for the correct number of input
	if(argc!=2) {
		Instructions();
		return 0;
	}

	srand(time(NULL));

	// Nuber of elements
	num_el = atoi(argv[1]);
	//printf("%d_%d_%d\n",num_el,rand(),RAND_MAX);

	//Creating password
	printf("Created password:\n\n");
	for(i=0; i<num_el; i++) {
		typ = rand() % ntype;
		r = rand() % tsize[typ] + initial[typ];
		printf("%c",r);
	}
	printf("\n\nby Marco Carmignotto, 07/02/2013\n\n");

	return 0;
}

void Instructions(void) {
	printf("\n\n\n");
	printf("**************************\n");
	printf("*  Password generator v1 *\n");
	printf("*  - Marco Carmignotto - *\n");
	printf("**************************\n");
	printf("\n\nUsage:\n");
	printf("   ./passgen <n>\n");
	printf("where:\n");
	printf("   <n> - number of characteres\n");
	printf("\n\n");

	return;
}
