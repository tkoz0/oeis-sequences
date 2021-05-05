/*
    Computes the difference of 2 sorted streams of nonnegative 64 bit integers.
    Output consists of "<%lu" lines for left stream and ">%lu" lines for right
    stream. The "<" or ">" indicates which stream it came from and %lu is the
    number.
    
    This program is slightly modified for the Fermat pseudoprimes use case. It
    expects both streams to end with a line containing "done".
*/

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{
    FILE *stream;
    uint64_t next;
    bool has_next;
}
STREAM_STATE;

// move stream state forward
static inline void advance(STREAM_STATE *ss)
{
    assert(ss->has_next);
    ss->has_next = fscanf(ss->stream,"%lu",&(ss->next)) == 1;
}

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr,"./a.out <left> <right>\n");
        fprintf(stderr,"use \"-\" for stdin\n");
        return 0;
    }
    if (!strcmp(argv[1],argv[2])) return 0; // same file
    STREAM_STATE L, R;
    L.stream = R.stream = NULL;
    if (!strcmp(argv[1],"-")) L.stream = stdin;
    else L.stream = fopen(argv[1],"r");
    if (!strcmp(argv[2],"-")) R.stream = stdin;
    else R.stream = fopen(argv[2],"r");
    if (!L.stream || !R.stream)
    {
        if (!L.stream) fprintf(stderr,"error opening: %s\n",argv[1]);
        if (!R.stream) fprintf(stderr,"error opening: %s\n",argv[2]);
        return 1;
    }
    // setup stream initial states
    L.has_next = fscanf(L.stream,"%lu",&L.next) == 1;
    R.has_next = fscanf(R.stream,"%lu",&R.next) == 1;
    while (L.has_next && R.has_next)
    {
        if (L.next < R.next)
        {
            printf("<%lu\n",L.next);
            advance(&L);
        }
        else if (R.next < L.next)
        {
            printf(">%lu\n",R.next);
            advance(&R);
        }
        else // equal, advance both
        {
            advance(&L);
            advance(&R);
        }
    }
    // finish unfinished stream
    while (L.has_next)
    {
        printf("<%lu\n",L.next);
        advance(&L);
    }
    while (R.has_next)
    {
        printf(">%lu\n",R.next);
        advance(&R);
    }
    char *doneL = NULL, *doneR = NULL;
    size_t lenL = 0, lenR = 0;
    // check for last line with "done"
    if (getline(&doneL,&lenL,L.stream) == -1
        || getline(&doneR,&lenR,R.stream) == -1)
        printf("error(reading)\n");
    else if (!strcmp(doneL,"done\n") && !strcmp(doneR,"done\n"))
        printf("done\n");
    else
        printf("error(values)\n");
    if (doneL) free(doneL);
    if (doneR) free(doneR);
    fclose(L.stream);
    fclose(R.stream);
    return 0;
}
