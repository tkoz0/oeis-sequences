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

Binary format for the recursion tree (-DWRITE_TREE):

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

Output format for statistics (-DWRITE_STATS):

    Comment lines start with #

    Begins with these 4 lines (value between < > are taken from arguments):
    # prime_type = <prime_type>
    # base = <base>
    # root = <root>
    # max_length = <max_length>

    May have this 5th line (for left-or-right truncatable primes):
    # NOTE: counts are not applicable

    Next is the CSV table, header is:
    digits,all,0,1,...,<max_children>
    <max_children> is determined from the base and prime type
    This splits the primes into based on how many child nodes they have

    For each nonzero digit length with at least 1 prime, 3 lines:
    - counts line (number of primes found)
    - min line (minimum of primes found)
    - max line (maximum of primes found)

    Finally, the hash line which can be used for verifying computation:
    # hash = <hash>

Hashing function (-DWRITE_STATS):

    Each node has a hash computed as follows (all calculations are 64 bit):
    node.value = the prime number of this node
    node.children = a list of child nodes

    rot32(n): return (n >> 32) | (n << 32)
    hash(node):
        h = (lower 64 bits of node.value) >> 1
        for child in node.children:
            d = path number for digit append
            c = hash(child)
            h ^= rot32(8191*(127*h - d) + c)
        return h

    The hash value output at the end is the hash of the root node
    
    The path numbers for digit appending (d above) are:
    - right/left truncatable: the digit
    - left-or-right truncatable:
      - the digit for left append
      - base + the digit for right append
    - left-and-right truncatable: (left digit)*base + (right digit)

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

/*
Ensure necessary requirements for preprocessor symbols
*/

#if defined (WRITE_TREE) && defined (WRITE_STATS)
#error "cannot define both WRITE_TREE and WRITE_STATS"
#endif

#if ! defined (WRITE_TREE) && ! defined (WRITE_STATS)
#error "must define either WRITE_TREE or WRITE_STATS"
#endif

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
#define BUFFER_SIZE (1<<20)
#ifdef WRITE_TREE
char *_g_buffer;
uint32_t _g_buffer_index;
#endif

// write to stdout for tree mode
static inline void write_buffer()
{
#ifdef WRITE_TREE
    if (write(1,_g_buffer,_g_buffer_index) != _g_buffer_index)
    {
        fprintf(stderr,"unable to write output\n");
        exit(1);
    }
    _g_buffer_index = 0;
#endif
}

// write a byte to the buffer, flush when full
static inline void write_byte(char b)
{
#ifdef WRITE_TREE
    _g_buffer[_g_buffer_index++] = b;
    if (_g_buffer_index == BUFFER_SIZE)
        write_buffer();
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

// macros for hashing the tree
#ifdef WRITE_TREE
typedef void tp_hash_t;
#endif
#ifdef WRITE_STATS
typedef uint64_t tp_hash_t;
#endif
#define HASH_INIT (mpz_get_ui(STACK_PREV)>>1)
#define HASH_DIGIT(h,d) (127*(h) - (d))
#define HASH_CHILD(h,c) (8191*(h) + (c))
#define HASH_SWAP(h) (((h) >> 32) | ((h) << 32))
#define HASH_UPDATE(h,d,c) ((h)^HASH_SWAP(HASH_CHILD(HASH_DIGIT(h,d),c)))

// right truncatable (A024770 for base 10)
tp_hash_t primes_r()
{
    ++_g_depth;
#ifdef WRITE_STATS
    uint32_t children = 0;
    tp_hash_t hash = HASH_INIT, child_hash;
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
                write_byte(d); // subtree
#ifdef WRITE_STATS
                child_hash =
#endif
                primes_r();
#ifdef WRITE_STATS
                ++children;
                hash = HASH_UPDATE(hash,d,child_hash);
#endif
            }
        }
    }
    --_g_depth;
    write_byte(255); // end
#ifdef WRITE_STATS
    ++_g_counts[_g_depth][children];
    update_min_max(children);
    return hash;
#endif
}

// left truncatable (A024785 for base 10)
tp_hash_t primes_l()
{
    ++_g_depth;
#ifdef WRITE_STATS
    uint32_t children = 0;
    tp_hash_t hash = HASH_INIT, child_hash;
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
                write_byte(d); // subtree
#ifdef WRITE_STATS
                child_hash =
#endif
                primes_l();
#ifdef WRITE_STATS
                ++children;
                hash = HASH_UPDATE(hash,d,child_hash);
#endif
            }
        }
    }
    --_g_depth;
    write_byte(255); // end
#ifdef WRITE_STATS
    ++_g_counts[_g_depth][children];
    update_min_max(children);
    return hash;
#endif
}

// left or right truncatable (A137812 for base 10)
tp_hash_t primes_lor()
{
    ++_g_depth;
#ifdef WRITE_STATS
    uint32_t children = 0;
    tp_hash_t hash = HASH_INIT, child_hash;
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
                write_byte(0); // subtree
                write_byte(d);
#ifdef WRITE_STATS
                child_hash =
#endif
                primes_lor();
#ifdef WRITE_STATS
                ++children;
                hash = HASH_UPDATE(hash,d,child_hash);
#endif
            }
        }
        // append right
        mpz_mul_ui(STACK_CURR,STACK_PREV,_g_base);
        for (uint32_t d = 1; d < _g_base; ++d)
        {
            mpz_add_ui(STACK_CURR,STACK_CURR,1);
            if (PRIME_TEST_CURR)
            {
                write_byte(1); // subtree
                write_byte(d);
#ifdef WRITE_STATS
                child_hash =
#endif
                primes_lor();
#ifdef WRITE_STATS
                ++children;
                hash = HASH_UPDATE(hash,_g_base+d,child_hash);
#endif
            }
        }
    }
    --_g_depth;
    write_byte(255); // end
#ifdef WRITE_STATS
    ++_g_counts[_g_depth][children];
    update_min_max(children);
    return hash;
#endif
}

// 2 digits at a time so use 2*_g_depth instead of _g_depth
#define LAR_POWER_INDEX (_g_rlen + (_g_depth << 1) - 1)

// left and right truncatable (A077390 for base 10)
tp_hash_t primes_lar()
{
    ++_g_depth;
#ifdef WRITE_STATS
    uint32_t children = 0;
    tp_hash_t hash = HASH_INIT, child_hash;
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
                    write_byte(dl); // subtree
                    write_byte(dr);
#ifdef WRITE_STATS
                    child_hash =
#endif
                    primes_lar();
#ifdef WRITE_STATS
                    ++children;
                    hash = HASH_UPDATE(hash,dl*_g_base+dr,child_hash);
#endif
                }
            }
            // backtrack right digit increment
            mpz_sub_ui(STACK_CURR,STACK_CURR,_g_base-1);
        }
    }
    --_g_depth;
    write_byte(255); // end
#ifdef WRITE_STATS
    ++_g_counts[_g_depth][children];
    update_min_max(children);
    return hash;
#endif
}

/*
Recursion setup functions
*/

// write the full tree given a root value
// set byte2 for cases where tree values are 2 bytes
// returns hash for the entire tree
tp_hash_t primes_init_root(tp_hash_t (*fptr)(), bool byte2)
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
#ifdef WRITE_STATS
    return
#endif
    fptr(); // writes subtrees and end
}

// write subtrees beginning with 1 digit
// specify -1 for byte2, otherwise specify the first byte value
// for left or right truncatable primes, it specifies side to append to
// for left and right truncatable primes, it should be 0
// updates the hash value with the subtrees and returns it
tp_hash_t primes_init_1digit(tp_hash_t (*fptr)(), int byte2, tp_hash_t *h0)
{
#ifdef WRITE_STATS
    tp_hash_t h = *h0, c;
#endif
    if (_g_maxlength < 1)
        return
#ifdef WRITE_STATS
        h
#endif
        ;
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
#ifdef WRITE_STATS
            c =
#endif
            fptr(); // subtrees and end
#ifdef WRITE_STATS
            h = HASH_UPDATE(h,root,c);
#endif
        }
    }
#ifdef WRITE_STATS
    return h;
#endif
}

// write subtrees beginning with 2 digits (2 bytes for root values)
// only used for left and right truncatable primes
// updates the hash value with the subtrees and returns it
tp_hash_t primes_init_2digit(tp_hash_t (*fptr)(), tp_hash_t *h0)
{
#ifdef WRITE_STATS
    tp_hash_t h = *h0, c;
#endif
    if (_g_maxlength < 2)
        return
#ifdef WRITE_STATS
        h
#endif
        ;
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
#ifdef WRITE_STATS
                c =
#endif
                fptr(); // subtrees and end
#ifdef WRITE_STATS
                h = HASH_UPDATE(h,rootl*_g_base+rootr,c);
#endif
            }
        }
#ifdef WRITE_STATS
    return h;
#endif
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
        if (count_all == 0) // skip rows with no primes
            continue;
        printf("%u,%lu",_g_rlen+i*mult,count_all);
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

static inline void write_hash(tp_hash_t *hash)
{
#ifdef WRITE_STATS
    printf("# hash = %lu\n",*hash);
#endif
}

/*
Recursion initialization functions
*/

#ifdef WRITE_TREE
#define HASH_ADDR NULL
#endif
#ifdef WRITE_STATS
#define HASH_ADDR (&hash)
#endif

// right truncatable (A024770 for base 10)
void primes_r_init()
{
#ifdef WRITE_STATS
    tp_hash_t hash = 0;
#endif
    if (mpz_cmp_ui(_g_root,0) > 0)
    {
//        if (!is_r_truncprime(_g_root,_g_base))
//        {
//            fprintf(stderr,"root is not right truncatable\n");
//            exit(1);
//        }
#ifdef WRITE_STATS
        hash =
#endif
        primes_init_root(primes_r,false);
    }
    else
    {
        write_byte(255); // root value
#ifdef WRITE_STATS
        hash =
#endif
        primes_init_1digit(primes_r,-1,HASH_ADDR);
        write_byte(255); // end
    }
    write_stats(1,true);
    write_hash(HASH_ADDR);
}

// left truncatable (A024785 for base 10)
void primes_l_init()
{
#ifdef WRITE_STATS
    tp_hash_t hash = 0;
#endif
    if (mpz_cmp_ui(_g_root,0) > 0)
    {
//        if (!is_l_truncprime(_g_root,_g_base))
//        {
//            fprintf(stderr,"root is not left truncatable\n");
//            exit(1);
//        }
#ifdef WRITE_STATS
        hash =
#endif
        primes_init_root(primes_l,false);
    }
    else
    {
        write_byte(255); // root value
#ifdef WRITE_STATS
        hash =
#endif
        primes_init_1digit(primes_l,-1,HASH_ADDR);
        write_byte(255); // end
    }
    write_stats(1,true);
    write_hash(HASH_ADDR);
}

// left or right truncatable (A137812 for base 10)
void primes_lor_init()
{
#ifdef WRITE_STATS
    tp_hash_t hash = 0;
#endif
    if (mpz_cmp_ui(_g_root,0) > 0)
    {
//        if (!is_lor_truncprime(_g_root,_g_base))
//        {
//            fprintf(stderr,"root is not left-or-right truncatable\n");
//            exit(1);
//        }
#ifdef WRITE_STATS
        hash =
#endif
        primes_init_root(primes_lor,true);
    }
    else
    {
        write_byte(255); // root value
        write_byte(255);
#ifdef WRITE_STATS
        hash =
#endif
        primes_init_1digit(primes_lor,0,HASH_ADDR);
        write_byte(255); // end
    }
    write_stats(1,true);
    write_hash(HASH_ADDR);
}

// left and right truncatable (A077390 for base 10)
void primes_lar_init()
{
#ifdef WRITE_STATS
    tp_hash_t hash = 0;
#endif
    if (mpz_cmp_ui(_g_root,0) > 0)
    {
//        if (!is_lar_truncprime(_g_root,_g_base))
//        {
//            fprintf(stderr,"root is not left-and-right truncatable\n");
//            exit(1);
//        }
#ifdef WRITE_STATS
        hash =
#endif
        primes_init_root(primes_lar,true);
        write_stats(2,true);
    }
    else
    {
        write_byte(255); // root value
        write_byte(255);
#ifdef WRITE_STATS
        hash =
#endif
        primes_init_1digit(primes_lar,0,HASH_ADDR);
        write_stats(2,true); // odd lengths
#ifdef WRITE_STATS
        clear_globals();
        init_globals();
        hash =
#endif
        primes_init_2digit(primes_lar,HASH_ADDR);
        write_stats(2,false); // even lengths
        write_byte(255); // end
    }
    write_hash(HASH_ADDR);
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
    write_buffer();
    mpz_clear(_g_root);
    clear_globals();
    return 0;
}
