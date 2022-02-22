/*
Generates truncatable primes.

-b base (--base base)
    number base to use for digit appends/truncations
    valid range is 2-255, default is 10
-l max_length (--max_length max_length)
    maximum length in digits to output, default is unlimited
-p type (--prime_type type)
    specify type of primes, currently supported are
    r - right truncatable (A024770 for base 10)
    l - left truncatable (A024785 for base 10)
    lor - left or right truncatable (A137812 for base 10)
    lar - left and right truncatable (A077390 for base 10)
    note that duplicates numbers may be included for "lor"
-r root (--root root)
    root number for recursion, as a 64 bit integer
    default is 0 (for recursing the entire tree)
    this option is intended for recursing subtrees on different threads
    it is not checked if root is a valid prime of the given type

Binary format for the recursion tree:

    supported up to base 255
    tree -> value [tree...] end
    value - the branch taken
        r,l - 1 byte with the digit appended (always nonzero)
        lor - 2 bytes, 0 for append left or 1 for right, followed by the digit
        lar - 2 bytes, the appended digits, left then right
              left digit is 0 for 1 digit roots
    tree... - zero or more trees
    end - a single 255 byte
    note that this format does not contain the base or root value
    the root value is 255 (1 or 2 bytes depending on prime type)
    
    pseudocode showing how to read this format
    next(bytes): extracts next byte from the stream
    peek(bytes): reads next byte from stream without extracting
    read_tree(bytes):
        value <- next(bytes) // or extract 2 bytes if appropriate
        while(peek(bytes) != 255)
            read_tree(bytes)
        end <- next(bytes)
        assert(end == 255)
*/

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

/*
Global constants and variables
*/

// command line arguments
const char *OPTION_STRING = "b:l:p:r:";
const struct option OPTION_LONG[] =
{
    { "base",       required_argument, 0, 'b' },
    { "max_length", required_argument, 0, 'l' },
    { "prime_type", required_argument, 0, 'p' },
    { "root",       required_argument, 0, 'r' },
    { 0,            0,                 0, 0   }
};
char *_g_prime_type;
mpz_t _g_root;

// buffer
#define BUFFER_SIZE (1<<16)
#ifdef WRITE_TREE
char *_g_buffer;
uint32_t _g_buffer_index;
#endif

// write a byte to the buffer, flush to stdout when full
static inline void write_byte(char b)
{
#ifdef WRITE_TREE
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
#endif
}

// for recursion
// _g_base must never be modified
// _g_depth must be set correctly before calling a recursive function
// _g_value must be returned to its original value afterward (backtracking)
// _g_powers must grow when necessary and _g_plen updated accordingly
uint32_t _g_base; // number base, used as a constant, must be >= 2
uint32_t _g_depth; // recursion depth
uint32_t _g_maxdepth; // depth limit
uint32_t _g_rlen; // root length (digits in specified base)
uint32_t _g_maxlength; // length limit (digits in specified base), constant
mpz_t *_g_stack; // recursion stack
uint32_t _g_slen; // length of _g_stack
mpz_t *_g_powers; // powers of the base, _g_powers[i] = base^i
uint32_t _g_plen; // length of _g_powers

// for stats output
// uses _g_slen to keep n digit min/max arrays the same length as _g_stack
// _g_pmin[i] = min prime of all on recursion level i
// _g_pmax[i] = max prime of all on recursion level i
// _g_counts[i][0] = count of non-leaf nodes on recursion level i
// _g_counts[i][1] = count of leaf nodes on recursion level i
// primes on recursion level i = _g_counts[i][0] + _g_counts[i][1]
#ifdef WRITE_STATS
mpz_t **_g_pmin, **_g_pmax;
uint64_t **_g_counts;
uint32_t _g_max_children;
#endif

// returns a pointer to the mpz_t variable representing base^p
// grows the array of base powers when necessary
static inline mpz_t *get_power(uint32_t p)
{
    if (p >= _g_plen)
    {
        _g_powers = realloc(_g_powers,sizeof(mpz_t)*(p+1));
        for (uint32_t i = _g_plen; i <= p; ++i)
        {
            mpz_init(_g_powers[i]);
            mpz_mul_ui(_g_powers[i],_g_powers[i-1],_g_base);
        }
        _g_plen = p+1;
    }
    return _g_powers+p;
}

// update the min and max using current stack value and recursion depth
// parameter is the number of child vertices
static inline void update_min_max(uint32_t cc)
{
#ifdef WRITE_STATS
    // also handles updating min for initial value 0
    if (mpz_cmp_ui(_g_pmin[_g_depth][cc],0) == 0
     || mpz_cmp(_g_pmin[_g_depth][cc],_g_stack[_g_depth]) > 0)
        mpz_set(_g_pmin[_g_depth][cc],_g_stack[_g_depth]);
    if (mpz_cmp(_g_pmax[_g_depth][cc],_g_stack[_g_depth]) < 0)
        mpz_set(_g_pmax[_g_depth][cc],_g_stack[_g_depth]);
#endif
}

// grows the stack to ensure space for _g_stack[i]
static inline void ensure_stack_space(uint32_t i)
{
    if (i >= _g_slen)
    {
        _g_stack = realloc(_g_stack,sizeof(*_g_stack)*(i+1));
#ifdef WRITE_STATS
        _g_pmin = realloc(_g_pmin,sizeof(*_g_pmin)*(i+1));
        _g_pmax = realloc(_g_pmax,sizeof(*_g_pmax)*(i+1));
        _g_counts = realloc(_g_counts,sizeof(*_g_counts)*(i+1));
#endif
        for (uint32_t j = _g_slen; j <= i; ++j)
        {
            mpz_init(_g_stack[j]);
#ifdef WRITE_STATS
            _g_pmin[j] = malloc(_g_max_children*sizeof(**_g_pmin));
            _g_pmax[j] = malloc(_g_max_children*sizeof(**_g_pmax));
            _g_counts[j] = calloc(_g_max_children,sizeof(**_g_counts));
            for (uint32_t k = 0; k < _g_max_children; ++k)
            {
                mpz_init(_g_pmin[j][k]);
                mpz_init(_g_pmax[j][k]);
            }
#endif
        }
        _g_slen = i+1;
    }
}

// checks if a string is a number
bool is_number(const char *s)
{
    if (!isdigit(*s)) // must start with digit
        return false;
    ++s;
    while (isdigit(*s)) // go past last digit
        ++s;
    return !(*s); // must be at null terminator
}

void init_globals()
{
#ifdef WRITE_TREE
    _g_buffer = malloc(BUFFER_SIZE);
    _g_buffer_index = 0;
#endif
#ifdef WRITE_STATS
    if (strcmp("lor",_g_prime_type) == 0)
        _g_max_children = 2*_g_base;
    else if (strcmp("lar",_g_prime_type) == 0)
        _g_max_children = _g_base*_g_base;
    else
        _g_max_children = _g_base;
    _g_pmin = malloc(sizeof(*_g_pmin));
    _g_pmax = malloc(sizeof(*_g_pmax));
    _g_counts = malloc(sizeof(*_g_counts));
    _g_pmin[0] = malloc(_g_max_children*sizeof(**_g_pmin));
    _g_pmax[0] = malloc(_g_max_children*sizeof(**_g_pmax));
    _g_counts[0] = calloc(_g_max_children,sizeof(**_g_counts));
    for (uint32_t k = 0; k < _g_max_children; ++k)
    {
        mpz_init(_g_pmin[0][k]);
        mpz_init(_g_pmax[0][k]);
    }
#endif
    _g_stack = malloc(sizeof(mpz_t));
    mpz_init(_g_stack[0]);
    _g_slen = 1;
    _g_powers = malloc(sizeof(mpz_t));
    mpz_init_set_ui(_g_powers[0],1);
    _g_plen = 1;
}

void clear_globals()
{
#ifdef WRITE_TREE
    free(_g_buffer);
#endif
    for (uint32_t i = 0; i < _g_slen; ++i)
    {
        mpz_clear(_g_stack[i]);
#ifdef WRITE_STATS
        for (uint32_t k = 0; k < _g_max_children; ++k)
        {
            mpz_clear(_g_pmin[i][k]);
            mpz_clear(_g_pmax[i][k]);
        }
        free(_g_pmin[i]);
        free(_g_pmax[i]);
        free(_g_counts[i]);
#endif
    }
#ifdef WRITE_STATS
    free(_g_pmin);
    free(_g_pmax);
    free(_g_counts);
#endif
    free(_g_stack);
    for (uint32_t i = 0; i < _g_plen; ++i)
        mpz_clear(_g_powers[i]);
    free(_g_powers);
}

/*
Truncatable primes recursion functions
These write the subtree bytes and the end byte
The caller is responsible for writing the value byte(s)
*/

// macros for recursion usage
#define STACK_CURR (_g_stack[_g_depth])
#define STACK_PREV (_g_stack[_g_depth-1])
#define POWER_CURR (*get_power(_g_rlen+_g_depth-1))
#define CHECK_STACK ensure_stack_space(_g_depth)

// primality test to use as a macro
// with GMP 6.2.0 it will run only a BPSW test when reps < 25
#define PRIME_TEST(n) mpz_probab_prime_p(n,0)
#define PRIME_TEST_CURR PRIME_TEST(STACK_CURR)

// right truncatable (A024770 for base 10)
void primes_r()
{
    ++_g_depth;
#ifdef WRITE_STATS
    uint32_t children = 0;
#endif
    if (_g_depth <= _g_maxdepth)
    {
        CHECK_STACK;
        // left shift to create a 0 digit on the right
        mpz_mul_ui(STACK_CURR,STACK_PREV,_g_base);
        for (uint32_t d = 1; d < _g_base; ++d)
        {
            // increment right digit
            mpz_add_ui(STACK_CURR,STACK_CURR,1);
            if (PRIME_TEST_CURR)
            {
#ifdef WRITE_STATS
                ++children;
#endif
                write_byte(d); // subtree
                primes_r();
            }
        }
    }
    --_g_depth;
#ifdef WRITE_STATS
    ++_g_counts[_g_depth][children];
    update_min_max(children);
#endif
    write_byte(255); // end
}

// left truncatable (A024785 for base 10)
void primes_l()
{
    ++_g_depth;
#ifdef WRITE_STATS
    uint32_t children = 0;
#endif
    if (_g_depth <= _g_maxdepth)
    {
        CHECK_STACK;
        mpz_set(STACK_CURR,STACK_PREV);
        for (uint32_t d = 1; d < _g_base; ++d)
        {
            // increment left digit
            mpz_add(STACK_CURR,STACK_CURR,POWER_CURR);
            if (PRIME_TEST_CURR)
            {
#ifdef WRITE_STATS
                ++children;
#endif
                write_byte(d); // subtree
                primes_l();
            }
        }
    }
    --_g_depth;
#ifdef WRITE_STATS
    ++_g_counts[_g_depth][children];
    update_min_max(children);
#endif
    write_byte(255); // end
}

// left or right truncatable (A137812 for base 10)
void primes_lor()
{
    ++_g_depth;
#ifdef WRITE_STATS
    uint32_t children = 0;
#endif
    if (_g_depth <= _g_maxdepth)
    {
        CHECK_STACK;
        // append left
        mpz_set(STACK_CURR,STACK_PREV);
        for (uint32_t d = 1; d < _g_base; ++d)
        {
            mpz_add(STACK_CURR,STACK_CURR,POWER_CURR);
            if (PRIME_TEST_CURR)
            {
#ifdef WRITE_STATS
                ++children;
#endif
                write_byte(0); // subtree
                write_byte(d);
                primes_lor();
            }
        }
        // append right
        mpz_mul_ui(STACK_CURR,STACK_PREV,_g_base);
        for (uint32_t d = 1; d < _g_base; ++d)
        {
            mpz_add_ui(STACK_CURR,STACK_CURR,1);
            if (PRIME_TEST_CURR)
            {
#ifdef WRITE_STATS
                ++children;
#endif
                write_byte(1); // subtree
                write_byte(d);
                primes_lor();
            }
        }
    }
    --_g_depth;
#ifdef WRITE_STATS
    ++_g_counts[_g_depth][children];
    update_min_max(children);
#endif
    write_byte(255); // end
}

// 2 digits at a time so use 2*_g_depth instead of _g_depth
#define LAR_POWER_INDEX (_g_rlen + (_g_depth << 1) - 1)

// left and right truncatable (A077390 for base 10)
void primes_lar()
{
    ++_g_depth;
#ifdef WRITE_STATS
    uint32_t children = 0;
#endif
    if ((_g_depth<<1) <= _g_maxdepth)
    {
        CHECK_STACK;
        mpz_mul_ui(STACK_CURR,STACK_PREV,_g_base); // shift left
        for (uint32_t dl = 1; dl < _g_base; ++dl)
        {
            // increment left digit
            mpz_add(STACK_CURR,STACK_CURR,*get_power(LAR_POWER_INDEX));
            // right digit loop
            for (uint32_t dr = 1; dr < _g_base; ++dr)
            {
                mpz_add_ui(STACK_CURR,STACK_CURR,1);
                if (PRIME_TEST_CURR)
                {
#ifdef WRITE_STATS
                    ++children;
#endif
                    write_byte(dl); // subtree
                    write_byte(dr);
                    primes_lar();
                }
            }
            // backtrack right digit increment
            mpz_sub_ui(STACK_CURR,STACK_CURR,_g_base-1);
        }
    }
    --_g_depth;
#ifdef WRITE_STATS
    ++_g_counts[_g_depth][children];
    update_min_max(children);
#endif
    write_byte(255); // end
}

/*
Recursion setup functions
*/

// write the full tree given a root value
// set byte2 for cases where tree values are 2 bytes
void primes_init_root(void (*fptr)(), bool byte2)
{
    write_byte(255); // root value
    if (byte2) // 2nd byte for root value
        write_byte(255);
    mpz_set(_g_stack[0],_g_root);
    _g_depth = 0;
    _g_rlen = 0;
    mpz_t root;
    mpz_init_set(root,_g_root);
    while (mpz_cmp_ui(root,0) > 0) // count digits
    {
        ++_g_rlen;
        mpz_div_ui(root,root,_g_base);
    }
    mpz_clear(root);
    _g_maxdepth = _g_maxlength - _g_rlen;
    fptr(); // writes subtrees and end
}

// write subtrees beginning with 1 digit
// specify -1 for byte2, otherwise specify the first byte value
// for left or right truncatable primes, it specifies side to append to
// for left and right truncatable primes, it should be 0
void primes_init_1digit(void (*fptr)(), int byte2)
{
    if (_g_maxlength < 1)
        return;
    for (uint64_t root = 2; root < _g_base; ++root)
    {
        mpz_set_ui(_g_stack[0],root);
        _g_depth = 0;
        _g_rlen = 1;
        _g_maxdepth = _g_maxlength - 1;
        if (PRIME_TEST(_g_stack[0]))
        {
            if (byte2 != -1)
                write_byte(byte2);
            write_byte(root); // root
            fptr(); // subtrees and end
        }
    }
}

// write subtrees beginning with 2 digits (2 bytes for root values)
// only used for left and right truncatable primes
void primes_init_2digit(void (*fptr)())
{
    if (_g_maxlength < 2)
        return;
    for (uint64_t rootl = 1; rootl < _g_base; ++rootl)
        for (uint64_t rootr = 0; rootr < _g_base; ++rootr)
        {
            mpz_set_ui(_g_stack[0],rootl*_g_base+rootr);
            _g_depth = 0;
            _g_rlen = 2;
            _g_maxdepth = _g_maxlength - 2;
            if (PRIME_TEST(_g_stack[0]))
            {
                write_byte(rootl); // root
                write_byte(rootr);
                fptr(); // subtrees and end
            }
        }
}

/*
Statistics output function
*/

// writes statistics in text format to stdout
// mult is number of digits appended at each recursion level
// header is whether to write first lines (not for 2nd time with lar primes)
void write_stats(uint32_t mult, bool header)
{
#ifdef WRITE_STATS
    if (header)
    {
        printf("# prime_type = %s\n",_g_prime_type);
        printf("# base = %u\n",_g_base);
        printf("# root = ");
        mpz_out_str(NULL,10,_g_root);
        printf("\n# max_length = %u\n",_g_maxlength);
        if (!strcmp(_g_prime_type,"lor"))
            printf("# NOTE: counts are not applicable\n");
        printf("digits,all");
        for (uint32_t k = 0; k < _g_max_children; ++k)
            printf(",%u",k);
        printf("\n");
    }
    uint64_t count_all;
    mpz_t min_all,max_all;
    mpz_init(min_all);
    mpz_init(max_all);
    for (uint32_t i = 0; i < _g_slen; ++i)
    {
        printf("%u,",_g_rlen+i*mult);
        count_all = 0;
        mpz_set_ui(min_all,0);
        mpz_set_ui(max_all,0);
        for (uint32_t k = 0; k < _g_max_children; ++k)
        {
            count_all += _g_counts[i][k];
            if (mpz_cmp_ui(min_all,0) == 0
             || (mpz_cmp(min_all,_g_pmin[i][k]) > 0
                 && mpz_cmp_ui(_g_pmin[i][k],0) != 0))
                mpz_set(min_all,_g_pmin[i][k]);
            if (mpz_cmp(max_all,_g_pmax[i][k]) < 0)
                mpz_set(max_all,_g_pmax[i][k]);
        }
        printf("%lu",count_all);
        for (uint32_t k = 0; k < _g_max_children; ++k)
            printf(",%lu",_g_counts[i][k]);
        printf("\n,");
        mpz_out_str(NULL,10,min_all);
        for (uint32_t k = 0; k < _g_max_children; ++k)
        {
            printf(",");
            mpz_out_str(NULL,10,_g_pmin[i][k]);
        }
        printf("\n,");
        mpz_out_str(NULL,10,max_all);
        for (uint32_t k = 0; k < _g_max_children; ++k)
        {
            printf(",");
            mpz_out_str(NULL,10,_g_pmax[i][k]);
        }
        printf("\n");
    }
    mpz_clear(min_all);
    mpz_clear(max_all);
#endif
}

/*
Recursion initialization functions
*/

// right truncatable (A024770 for base 10)
void primes_r_init()
{
    if (mpz_cmp_ui(_g_root,0) > 0)
    {
        if (!is_r_truncprime(_g_root,_g_base))
        {
            fprintf(stderr,"root is not right truncatable\n");
            exit(0);
        }
        primes_init_root(primes_r,false);
    }
    else
    {
        write_byte(255); // root value
        primes_init_1digit(primes_r,-1);
        write_byte(255); // end
    }
    write_stats(1,true);
}

// left truncatable (A024785 for base 10)
void primes_l_init()
{
    if (mpz_cmp_ui(_g_root,0) > 0)
    {
        if (!is_l_truncprime(_g_root,_g_base))
        {
            fprintf(stderr,"root is not left truncatable\n");
            exit(0);
        }
        primes_init_root(primes_l,false);
    }
    else
    {
        write_byte(255); // root value
        primes_init_1digit(primes_l,-1);
        write_byte(255); // end
    }
    write_stats(1,true);
}

// left or right truncatable (A137812 for base 10)
void primes_lor_init()
{
    if (mpz_cmp_ui(_g_root,0) > 0)
    {
        if (!is_lor_truncprime(_g_root,_g_base))
        {
            fprintf(stderr,"root is not left-or-right truncatable\n");
            exit(0);
        }
        primes_init_root(primes_lor,true);
    }
    else
    {
        write_byte(255); // root value
        write_byte(255);
        primes_init_1digit(primes_lor,0);
        write_byte(255); // end
    }
    write_stats(1,true);
}

// left and right truncatable (A077390 for base 10)
void primes_lar_init()
{
    if (mpz_cmp_ui(_g_root,0) > 0)
    {
        if (!is_lar_truncprime(_g_root,_g_base))
        {
            fprintf(stderr,"root is not left-and-right truncatable\n");
            exit(0);
        }
        primes_init_root(primes_lar,true);
    }
    else
    {
        write_byte(255); // root value
        write_byte(255);
        primes_init_1digit(primes_lar,0);
        write_stats(2,true); // odd lengths
#ifdef WRITE_STATS
        clear_globals();
        init_globals();
#endif
        primes_init_2digit(primes_lar);
        write_stats(2,false); // even lengths
        write_byte(255); // end
    }
}

/*
Main function
*/

int main(int argc, char **argv)
{
    // set default values
    _g_base = 10;
    _g_maxlength = -1;
    _g_prime_type = NULL;
    mpz_init(_g_root);
    if (argc < 2)
    {
        fprintf(stderr,"truncprimes <-p prime_type> "
                        "[-b base] [-l max_length] [-r root]\n");
        return 0;
    }
    // read options
    int o;
    while ((o = getopt_long(argc,argv,OPTION_STRING,OPTION_LONG,NULL)) != -1)
    {
        switch (o)
        {
        case 'b': // base
            if (!is_number(optarg))
            {
                fprintf(stderr,"base must be a number\n");
                return 0;
            }
            _g_base = atoi(optarg);
            break;
        case 'l': // max length
            if (!is_number(optarg))
            {
                fprintf(stderr,"max length must be a number\n");
                return 0;
            }
            _g_maxlength = atoi(optarg);
            break;
        case 'p': // prime type
            _g_prime_type = optarg;
            break;
        case 'r': // root
            if (mpz_set_str(_g_root,optarg,10) != 0
             || mpz_cmp_ui(_g_root,0) < 0)
            {
                fprintf(stderr,"root must be a nonnegative integer\n");
                return 0;
            }
            break;
        default:
            fprintf(stderr,"error parsing arguments\n");
            fprintf(stderr,"truncprimes <-p prime_type> "
                        "[-b base] [-l max_length] [-r root]\n");
            return 0;
        }
    }
    if (_g_base < 2 || _g_base > 255)
    {
        fprintf(stderr,"base %u out of valid range (2-255)\n",_g_base);
        return 0;
    }
    if (!_g_prime_type)
    {
        fprintf(stderr,"must specify prime type\n");
        return 0;
    }
    init_globals();
    if (strcmp(_g_prime_type,"r") == 0)
        primes_r_init();
    else if (strcmp(_g_prime_type,"l") == 0)
        primes_l_init();
    else if (strcmp(_g_prime_type,"lor") == 0)
        primes_lor_init();
    else if (strcmp(_g_prime_type,"lar") == 0)
        primes_lar_init();
    else
        fprintf(stderr,"invalid prime type: %s\n",_g_prime_type);
    // flush buffer and exit
#ifdef WRITE_TREE
    if (write(1,_g_buffer,_g_buffer_index) != _g_buffer_index)
    {
        fprintf(stderr,"unable to write output\n");
        exit(1);
    }
#endif
    mpz_clear(_g_root);
    clear_globals();
    return 0;
}
