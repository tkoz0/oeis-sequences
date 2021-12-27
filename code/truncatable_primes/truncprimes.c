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
    value - the branch taken, or zero byte(s) for the root
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

// buffer
#define BUFFER_SIZE (1<<16)
char *_g_buffer;
uint32_t _g_buffer_index;

// write a byte to the buffer, flush to stdout when full
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

// grows the stack to ensure space for _g_stack[i]
static inline void ensure_stack_space(uint32_t i)
{
    if (i >= _g_slen)
    {
        _g_stack = realloc(_g_stack,sizeof(mpz_t)*(i+1));
        for (uint32_t j = _g_slen; j <= i; ++j)
            mpz_init(_g_stack[j]);
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
    _g_buffer = malloc(BUFFER_SIZE);
    _g_buffer_index = 0;
    _g_stack = malloc(sizeof(mpz_t));
    mpz_init(_g_stack[0]);
    _g_slen = 1;
    _g_powers = malloc(sizeof(mpz_t));
    mpz_init_set_ui(_g_powers[0],1);
    _g_plen = 1;
}

void clear_globals()
{
    free(_g_buffer);
    for (uint32_t i = 0; i < _g_slen; ++i)
        mpz_clear(_g_stack[i]);
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
// with GMP 6.2.0 it will run a BPSW test only if reps < 25
#define PRIME_TEST(n) mpz_probab_prime_p(n,0)
#define PRIME_TEST_CURR PRIME_TEST(STACK_CURR)

// right truncatable (A024770 for base 10)
void primes_r()
{
    ++_g_depth;
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
                primes_r();
            }
        }
    }
    --_g_depth;
    write_byte(255); // end
}

// left truncatable (A024785 for base 10)
void primes_l()
{
    ++_g_depth;
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
                primes_l();
            }
        }
    }
    --_g_depth;
    write_byte(255); // end
}

// left or right truncatable (A137812 for base 10)
void primes_lor()
{
    ++_g_depth;
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
                write_byte(1); // subtree
                write_byte(d);
                primes_lor();
            }
        }
    }
    --_g_depth;
    write_byte(255); // end
}

// 2 digits at a time so use 2*_g_depth instead of _g_depth
#define LAR_POWER_INDEX (_g_rlen + (_g_depth << 1) - 1)

// left and right truncatable (A077390 for base 10)
void primes_lar()
{
    ++_g_depth;
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
                    primes_lar();
                }
            }
            // backtrack right digit increment
            mpz_sub_ui(STACK_CURR,STACK_CURR,_g_base-1);
        }
    }
    --_g_depth;
    write_byte(255); // end
}

/*
Recursion setup functions
*/

// write the full tree given a root value
// set byte2 for cases where tree values are 2 bytes
void primes_init_root(uint64_t root, void (*fptr)(), bool byte2)
{
    write_byte(255); // root value
    if (byte2) // 2nd byte for root value
        write_byte(255);
    mpz_set_ui(_g_stack[0],root);
    _g_depth = 0;
    _g_rlen = 0;
    _g_maxdepth = _g_maxlength;
    while (root)
        ++_g_rlen, root /= _g_base;
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

// right truncatable (A024770 for base 10)
void primes_r_init(uint64_t root)
{
    if (root)
        primes_init_root(root,primes_r,false);
    else
    {
        write_byte(255); // root value
        primes_init_1digit(primes_r,-1);
        write_byte(255); // end
    }
}

// left truncatable (A024785 for base 10)
void primes_l_init(uint64_t root)
{
    if (root)
        primes_init_root(root,primes_l,false);
    else
    {
        write_byte(255); // root value
        primes_init_1digit(primes_l,-1);
        write_byte(255); // end
    }
}

// left or right truncatable (A137812 for base 10)
void primes_lor_init(uint64_t root)
{
    if (root)
        primes_init_root(root,primes_lor,true);
    else
    {
        write_byte(255); // root value
        write_byte(255);
        primes_init_1digit(primes_lor,0);
        write_byte(255); // end
    }
}

// left and right truncatable (A077390 for base 10)
void primes_lar_init(uint64_t root)
{
    if (root)
        primes_init_root(root,primes_lar,true);
    else
    {
        write_byte(255); // root value
        write_byte(255);
        primes_init_1digit(primes_lar,0);
        primes_init_2digit(primes_lar);
        write_byte(255); // end
    }
}

/*
Main function
*/

int main(int argc, char **argv)
{
    init_globals();
    // set default values
    _g_base = 10;
    _g_maxlength = -1;
    char *prime_type = NULL;
    uint64_t root = 0;
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
            prime_type = optarg;
            break;
        case 'r': // root
            if (!is_number(optarg))
            {
                fprintf(stderr,"root must be a number\n");
                return 0;
            }
            root = atoll(optarg);
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
    if (!prime_type)
    {
        fprintf(stderr,"must specify prime type\n");
        return 0;
    }
    if (strcmp(prime_type,"r") == 0)
        primes_r_init(root);
    else if (strcmp(prime_type,"l") == 0)
        primes_l_init(root);
    else if (strcmp(prime_type,"lor") == 0)
        primes_lor_init(root);
    else if (strcmp(prime_type,"lar") == 0)
        primes_lar_init(root);
    else
        fprintf(stderr,"invalid prime type: %s\n",prime_type);
    // flush buffer and exit
    if (write(1,_g_buffer,_g_buffer_index) != _g_buffer_index)
    {
        fprintf(stderr,"unable to write output\n");
        exit(1);
    }
    clear_globals();
    return 0;
}
