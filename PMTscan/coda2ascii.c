/************************************************* 
*    Simple test program for opening and reading *
*   an existing CODA event format file.          *
*                                                *
*   Author: David Abbott  CEBAF                  *
*                                                *
*   Arguments to the routine are:                *
*                                                *
*             evt <filename>                     *
*                                                *
**************************************************/

#include <stdio.h>

#define MAXBUFLEN  32768

main (argc, argv)
     int argc;
     char **argv;
{
  int handle, physEvtNum, status;
  unsigned long evtype;
  unsigned long buf[MAXBUFLEN];
  
  if(argc != 2) {
    printf("Incorrect number of arguments:\n");
    printf("  usage: evt filename\n");
    exit(-1);
  }

  if ( (status = evOpen(argv[1],"r",&handle)) < 0) {
    printf("Unable to open file %s status = 0x%x\n",argv[1],status);
    exit(-1);
  }

  while ((status=evRead(handle,buf,MAXBUFLEN))==0) {
    evtype = (buf[1]&0xffff0000)>>16;
    if (evtype < 16) { /* Physics Event */
      int dataStart = 8, dataEnd = buf[7]+7;
      while (dataStart <= dataEnd){
	unsigned long dataWord = buf[dataStart++];
	
	printf("%x,", dataWord);
      }
      printf("\n");
    }
  }

  evClose(handle);
  
  exit(0);
}
