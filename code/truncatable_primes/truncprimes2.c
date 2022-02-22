#include <assert.h>
#include <ctype.h>
#include <getopt.h>
#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "tp_util.h"

// command line
const char *OPTION_STRING = "b:l:o:p:r:";
const struct option OPTION_LONG[] =
{
    {"base",required_argument,0,'b'},
    {"maxlength",required_argument,0,'l'},
    {"output",required_argument,0,'o'},
    {"primetype",required_argument,0,'p'},
    {"root",required_argument,0,'r'},
    {0,0,0,0}
};

#define BUFFER_SIZE (1<<20)
char *_g_buffer;
uint32_t _g_buffer_index;

static inline void write_byte(char b)
{
    _g_buffer[_g_buffer_index++] = b;
    if (_g_buffer_index == BUFFER_SIZE)
    {
        if (write(1,_g_buffer,BUFFER_SIZE) != BUFFER_SIZE)
        {
            fprintf(stderr,"unable to write output\n");
            exit(1);
        }
        _g_buffer_index = 0;
    }
}

static inline void write_bytes(const char *v, uint32_t l)
{
    for (uint32_t i = 0; i < l; ++i)
        write_byte(v[i]);
}

static inline void flush_bytes()
{
    if (write(1,_g_buffer,_g_buffer_index) != _g_buffer_index)
    {
        fprintf(stderr,"unable to write output\n");
        exit(1);
    }
    _g_buffer_index = 0;
}

// global variables for arguments and stuff
uint32_t arg_base;
uint32_t arg_maxlength;
char *arg_primetype;
char *arg_output;
mpz_t arg_root;
bool stats_output;
uint32_t (*TP_next)(TP_STATE*,char*,TP_VALUE*);

// global variables for stats
// index with length, then number of children
mpz_t **s_min, **s_max;
uint64_t **s_count;
uint32_t s_length, s_maxchildren;

void resize_stats(uint32_t length)
{
    assert(length > s_length);
    if (s_length == 0)
    {
        s_min = malloc(length*sizeof(*s_min));
        s_max = malloc(length*sizeof(*s_max));
        s_count = malloc(length*sizeof(*s_count));
    }
    else
    {
        s_min = realloc(s_min,length*sizeof(*s_min));
        s_max = realloc(s_max,length*sizeof(*s_max));
        s_count = realloc(s_count,length*sizeof(*s_count));
    }
    while (s_length < length)
    {
        s_min[s_length] = malloc(s_maxchildren*sizeof(**s_min));
        s_max[s_length] = malloc(s_maxchildren*sizeof(**s_max));
        for (uint32_t i = 0; i < s_maxchildren; ++i)
        {
            mpz_init(s_min[s_length][i]);
            mpz_init(s_max[s_length][i]);
        }
        s_count[s_length] = calloc(s_maxchildren,sizeof(**s_count));
        ++s_length;
    }
}

void clear_stats()
{
    if (!s_length)
        return;
    for (uint32_t i = 0; i < s_length; ++i)
    {
        for (uint32_t j = 0; j < s_maxchildren; ++j)
        {
            mpz_clear(s_min[i][j]);
            mpz_clear(s_max[i][j]);
        }
        free(s_min[i]);
        free(s_max[i]);
        free(s_count[i]);
    }
    free(s_min);
    free(s_max);
    free(s_count);
    s_length = 0;
}

// for bytes output, writes the tree for the given root
// for stats output, updates global stats objects
void process_root(const mpz_t root, const char *rootbytes)
{
    TP_STATE state;
    TP_VALUE value;
    char bytes[2];
    uint32_t ret;
    TP_init(&state,arg_base,root,arg_maxlength,rootbytes,
        stats_output ? TP_POST_ORDER : TP_BYTES_ONLY);
    while ((ret = TP_next(&state,bytes,&value)))
    {
        if (stats_output)
        {
            if (value.len)
            {
                ;
            }
// TODO update stats
            fprintf(stderr,"stats not implemented\n");
            exit(1);
        }
        else
            write_bytes(bytes,ret);
    }
    TP_clear(&state);
}

int main(int argc, char **argv)
{
    _g_buffer = malloc(BUFFER_SIZE);
    _g_buffer_index = 0;
    // default values
    arg_base = 10; // base
    arg_maxlength = -1; // recursion depth limit
    arg_primetype = NULL; // "r","l","lor","lar"
    arg_output = NULL; // "tree","stats"
    mpz_init(arg_root);
    // help message
    if (argc < 2)
    {
        fprintf(stderr,"help\n");
        return 0;
    }
    // read options
    int o;
    while ((o = getopt_long(argc,argv,OPTION_STRING,OPTION_LONG,NULL)) != -1)
    {
        switch (o)
        {
        case 'b': // base
            arg_base = atoi(optarg);
            break;
        case 'l': // maxlength
            arg_maxlength = atoi(optarg);
            break;
        case 'o': // output
            arg_output = optarg;
            break;
        case 'p': // primetype
            arg_primetype = optarg;
            break;
        case 'r': // root
            if (mpz_set_str(arg_root,optarg,10) != 0
             || mpz_cmp_ui(arg_root,0) < 0)
            {
                fprintf(stderr,"root must be integer >= 0\n");
                return 0;
            }
            break;
        default:
            fprintf(stderr,"error parsing arguments\n");
            return 0;
        }
    }
    if (arg_base < 2 || arg_base > 255)
    {
        fprintf(stderr,"base (%u) must be in range 2-255\n",arg_base);
        return 0;
    }
    if (!arg_primetype)
    {
        fprintf(stderr,"must specify prime type\n");
        return 0;
    }
    if (!arg_output)
    {
        fprintf(stderr,"must specify output type\n");
        return 0;
    }
    // set function to use
    if (!strcmp(arg_primetype,"r"))
        TP_next = TP_next_r;
    else if (!strcmp(arg_primetype,"l"))
        TP_next = TP_next_l;
    else if (!strcmp(arg_primetype,"lor"))
        TP_next = TP_next_lor;
    else if (!strcmp(arg_primetype,"lar"))
        TP_next = TP_next_lar;
    else
    {
        fprintf(stderr,"invalid prime type %s\n",arg_primetype);
        return 0;
    }
    // determine output type
    if (!strcmp(arg_output,"tree"))
        stats_output = false;
    else if (!strcmp(arg_output,"stats"))
        stats_output = true;
    else
    {
        fprintf(stderr,"invalid output type %s\n",arg_output);
        return 0;
    }
    s_maxchildren = arg_base;
    if (!strcmp(arg_primetype,"lor"))
        s_maxchildren *= 2;
    else if (!strcmp(arg_primetype,"lar"))
        s_maxchildren *= arg_base;
    s_length = 0;
    // write results
    char rootbytes[2];
    if (mpz_cmp_ui(arg_root,0) == 0)
    {
        uint32_t limit = arg_base;
        mpz_t root;
        mpz_init(root);
        if (!strcmp(arg_primetype,"lar"))
            limit *= arg_base;
        write_byte(255); // root \xff or \xff\xff
        if (!strcmp(arg_primetype,"lor") || !strcmp(arg_primetype,"lar"))
            write_byte(255);
        for (uint32_t r = 1; r < limit; ++r)
        {
            mpz_set_ui(root,r);
            if (!PRIME_TEST(root))
                continue;
            if (!strcmp(arg_primetype,"lor"))
                rootbytes[0] = 0, rootbytes[1] = r;
            else if (!strcmp(arg_primetype,"lar"))
                rootbytes[0] = r/arg_base, rootbytes[1] = r%arg_base;
            else
                rootbytes[0] = r;
            process_root(root,rootbytes);
        }
        write_byte(255); // end
        mpz_clear(root);
    }
    else // single subtree only
    {
        rootbytes[0] = rootbytes[1] = 255; // \xff\xff root value
        process_root(arg_root,rootbytes);
    }
    if (stats_output)
    {
// TODO write stats
        fprintf(stderr,"stats not implemented\n");
        return 0;
    }
    else
        flush_bytes();
    mpz_clear(arg_root);
    free(_g_buffer);
    return 0;
}

// TODO change root value to a special exception with prime type and root
// store this information in ascii followed by some null terminator maybe
// then recursively read child subtrees and end
