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
        r,l - 1 byte with the digit appended
        lor - 2 bytes, 1 for append right or 2 for left, followed by the digit
        lar - 2 bytes, the appended digits, left then right
    tree... - zero or more trees
    end - a single 0 byte
    note that this format does not contain the base or root value
    to avoid conflict with end byte, all values must start with a nonzero byte
    (only the root value is permitted to start with a 0 byte)
    
    pseudocode showing how to read this format
    next(bytes): extracts next byte from the stream
    peek(bytes): reads next byte from stream without extracting
    read_tree(bytes):
        value <- next(bytes) // or extract 2 bytes if appropriate
        while(peek(bytes) != 0)
            read_tree(bytes)
        end <- next(bytes)
        assert(end == 0)
*/

#include <assert.h>
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
#define BUFFER_SIZE (1<<20)
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

// integers 0-3 for where basic C type is not supported (such as exponent base)
mpz_t MPZ_0, MPZ_1, MPZ_2, MPZ_3;

// for temporary variables in functions
mpz_t _g_tmp0, _g_tmp1;

// for recursion
// _g_base must never be modified
// _g_depth must be set correctly before calling a recursive function
// _g_value must be returned to its original value afterward (backtracking)
// _g_powers must grow when necessary and _g_plen updated accordingly
uint32_t _g_base; // number base, used as a constant, must be >= 2
uint32_t _g_depth; // recursion depth / number of digits (in given base)
uint32_t _g_maxdepth; // depth limit
mpz_t _g_value; // number used in recursion, manipulated in place
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

void init_globals()
{
    _g_buffer = malloc(BUFFER_SIZE);
    _g_buffer_index = 0;
    mpz_init_set_ui(MPZ_0,0);
    mpz_init_set_ui(MPZ_1,1);
    mpz_init_set_ui(MPZ_2,2);
    mpz_init_set_ui(MPZ_3,3);
    mpz_inits(_g_tmp0,_g_tmp1,NULL);
    mpz_init(_g_value);
    _g_powers = malloc(sizeof(mpz_t));
    mpz_init_set_ui(_g_powers[0],1);
    _g_plen = 1;
}

void clear_globals()
{
    free(_g_buffer);
    mpz_clears(MPZ_0,MPZ_1,MPZ_2,MPZ_3,NULL);
    mpz_clears(_g_tmp0,_g_tmp1,NULL);
    mpz_clear(_g_value);
    for (uint32_t i = 0; i < _g_plen; ++i)
        mpz_clear(_g_powers[i]);
    free(_g_powers);
}

/*
Primality testing
*/

// PRP(2) test
// should avoid calling with even inputs for efficiency
// must not call with n == 2
static inline bool is_prime_prp2(const mpz_t n)
{
    assert(mpz_odd_p(n));
    mpz_sub_ui(_g_tmp0,n,1); // n-1
    mpz_powm(_g_tmp0,MPZ_2,_g_tmp0,n); // 2^(n-1) mod n
    return mpz_cmp_ui(_g_tmp0,1) == 0; // == 1
}

// SPRP(2) test
// must only call with odd n > 2
static inline bool is_prime_sprp2(const mpz_t n)
{
    assert(mpz_odd_p(n));
    mpz_sub_ui(_g_tmp0,n,1); // _g_tmp0 == n-1
    uint32_t s = mpz_scan1(_g_tmp0,0);
    assert(s);
    mpz_fdiv_q_2exp(_g_tmp1,_g_tmp0,s); // n-1 == _g_tmp1 * 2^s
    mpz_powm(_g_tmp1,MPZ_2,_g_tmp1,n);
    if (mpz_cmp_ui(_g_tmp1,1) == 0 || mpz_cmp(_g_tmp1,_g_tmp0) == 0) // 1,n-1
        return true;
    --s;
    while (s--) // s-1 squarings
    {
        mpz_powm_ui(_g_tmp1,_g_tmp1,2,n); // square
        if (mpz_cmp(_g_tmp1,_g_tmp0) == 0) // n-1
            return true;
        // squaring 0 or 1 further will never result in n-1
        if (mpz_cmp_ui(_g_tmp1,1) <= 0) // from GMP code
            return false;
    }
    return false;
}

// BPSW test (not implemented)
static inline bool is_prime_bpsw(const mpz_t n)
{
    assert(0);
    return false;
}

/*
Truncatable primes recursion functions
These write the subtree bytes and the end byte
The caller is responsible for writing the value byte(s)
*/

// primality test to use as a macro, n=2 must be handled separately
#define PRIME_TEST(n) (mpz_odd_p(n) && is_prime_sprp2(n))

void primes_r()
{
    ++_g_depth;
    if (_g_depth <= _g_maxdepth)
    {
        // left shift to create a 0 digit on the right
        mpz_mul_ui(_g_value,_g_value,_g_base);
        for (uint32_t d = 1; d < _g_base; ++d)
        {
            // increment right digit
            mpz_add_ui(_g_value,_g_value,1);
            if (PRIME_TEST(_g_value))
            {
                write_byte(d); // subtree
                primes_r();
            }
        }
        // backtrack
        mpz_div_ui(_g_value,_g_value,_g_base);
    }
    --_g_depth;
    write_byte(0); // end
}

void primes_l()
{
    ++_g_depth;
    if (_g_depth <= _g_maxdepth)
    {
        for (uint32_t d = 1; d < _g_base; ++d)
        {
            // increment left digit
            mpz_add(_g_value,_g_value,*get_power(_g_depth-1));
            if (PRIME_TEST(_g_value))
            {
                write_byte(d); // subtree
                primes_l();
            }
        }
        // backtrack
        mpz_submul_ui(_g_value,*get_power(_g_depth-1),_g_base-1);
    }
    --_g_depth;
    write_byte(0); // end
}

void primes_lor()
{
    ;
}

void primes_lar()
{
    ;
}

/*
Recursion setup functions
*/

void primes_r_init(uint64_t root)
{
    write_byte(0); // root value
    if (root)
    {
        mpz_set_ui(_g_value,root);
        _g_depth = 0;
        while (root) // count digits
            ++_g_depth, root /= _g_base;
        primes_r();
        return;
    }
    // root is 0, initialize with single digit primes
    if (1 <= _g_maxdepth)
        for (root = 2; root < _g_base; ++root)
        {
            mpz_set_ui(_g_value,root);
            _g_depth = 1;
            if (mpz_cmp_ui(_g_value,2) == 0 || PRIME_TEST(_g_value))
            {
                write_byte(root); // subtree
                primes_r();
            }
        }
    write_byte(0); // end
}

void primes_l_init(uint64_t root)
{
    write_byte(0); // root value
    if (root)
    {
        mpz_set_ui(_g_value,root);
        _g_depth = 0;
        while (root) // count digits
            ++_g_depth, root /= _g_base;
        primes_l();
        return;
    }
    // root is 0, initialize with single digit primes
    if (1 <= _g_maxdepth)
        for (root = 2; root < _g_base; ++root)
        {
            mpz_set_ui(_g_value,root);
            _g_depth = 1;
            if (mpz_cmp_ui(_g_value,2) == 0 || PRIME_TEST(_g_value))
            {
                write_byte(root); // subtree
                primes_l();
            }
        }
    write_byte(0); // end
}

void primes_lor_init(uint64_t root)
{
    ;
}

void primes_lar_init(uint64_t root)
{
    ;
}

/*
Main function
*/

int main(int argc, char **argv)
{
    init_globals();
    // set default values
    _g_base = 10;
    _g_maxdepth = -1;
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
            _g_base = atoi(optarg);
            break;
        case 'l': // max length
            _g_maxdepth = atoi(optarg);
            break;
        case 'p': // prime type
            prime_type = optarg;
            break;
        case 'r': // root
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
    {
        fprintf(stderr,"not implemented\n");
        primes_lor_init(root);
    }
    else if (strcmp(prime_type,"lar") == 0)
    {
        fprintf(stderr,"not implemented\n");
        primes_lar_init(root);
    }
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
