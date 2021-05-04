/*
    Computes the difference of 2 sorted streams of nonnegative 64 bit integers.
    Output consists of "<%lu" lines for left stream and ">%lu" lines for right
    stream. The "<" or ">" indicates which stream it came from and %lu is the
    number. This does runtime assertion of sorted order.
*/

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr,"./a.out <left> <right>\n");
        fprintf(stderr,"use \"-\" for stdin\n");
        return 0;
    }
    if (!strcmp(argv[1],argv[2])) return 0; // same file
    FILE *L, *R;
    if (!strcmp(argv[1],"-")) L = stdin;
    else L = fopen(argv[1],"r");
    if (!strcmp(argv[2],"-")) R = stdin;
    else R = fopen(argv[2],"r");
    if (!L || !R)
    {
        fprintf(stderr,"cannot open both files\n");
        return 1;
    }
    uint64_t l = 0, r = 0, tmpl, tmpr; // values in the streams
    // read next integers for each iteration
    while (fscanf(L,"%lu",&tmpl) == 1 && fscanf(R,"%lu",&tmpr) == 1)
    {
        assert(tmpl >= l);
        assert(tmpr >= r);
        l = tmpl;
        r = tmpr;
        // repeat until equal, then outer loop reads the next ones
        while (l != r)
        {
            while (l < r)
            {
                printf("<%lu\n",l);
                if (fscanf(L,"%lu",&tmpl) != 1) goto eof;
                assert(tmpl >= l);
                l = tmpl;
            }
            while (r < l)
            {
                printf(">%lu\n",r);
                if (fscanf(R,"%lu",&tmpr) != 1) goto eof;
                assert(tmpr >= r);
                r = tmpr;
            }
        }
    }
    eof: // go here once a file reaches EOF
    // complete the stream not at eof
    while (fscanf(L,"%lu",&tmpl) == 1)
    {
        assert(tmpl >= l);
        l = tmpl;
        printf("<%lu\n",l);
    }
    while (fscanf(R,"%lu",&tmpr) == 1)
    {
        assert(tmpr >= r);
        r = tmpr;
        printf(">%lu\n",r);
    }
    return 0;
}
