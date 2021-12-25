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
    the root value is 0 (or 2 zero bytes)
    
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

// Modulus for small primes test, need to be able to multiply digit (up to 254)
// Moduli for 32 and 64 bit respectively
#define PRIMORIAL_19 9699690u
#define PRIMORIAL_43 13082761331670030uL
// The small primes modulus to use
#define SPMOD PRIMORIAL_43
// Type to use for small primes division
typedef uint64_t spmod_t;
// Limit for proving primality with division by these primes
#define TDIV_LIMIT (47*47)

// integers 0-3 for where basic C type is not supported (such as exponent base)
mpz_t MPZ_0, MPZ_1, MPZ_2, MPZ_3;

// for recursion
// _g_base must never be modified
// _g_depth must be set correctly before calling a recursive function
// _g_value must be returned to its original value afterward (backtracking)
// _g_powers must grow when necessary and _g_plen updated accordingly
uint32_t _g_base; // number base, used as a constant, must be >= 2
uint32_t _g_depth; // recursion depth
uint32_t _g_maxdepth; // recursion depth limit
uint32_t _g_rlen; // root length (digits in specified base)
uint32_t _g_maxlength; // length limit (digits in specified base), constant
mpz_t *_g_stack; // recursion stack
spmod_t *_g_spmods; // other component of recursion stack
uint32_t _g_slen; // length of _g_stack
mpz_t *_g_powers; // powers of the base, _g_powers[i] = base^i
spmod_t *_g_power_spmods; // powers of base modulo SPMOD
uint32_t _g_plen; // length of _g_powers and _g_power_spmods

// returns a pointer to the mpz_t variable representing base^p
// grows the array of base powers when necessary
// also grows the _g_power_spmods array
static inline mpz_t *get_power(uint32_t p)
{
    if (p >= _g_plen)
    {
        _g_powers = realloc(_g_powers,sizeof(mpz_t)*(p+1));
        _g_power_spmods = realloc(_g_power_spmods,sizeof(spmod_t)*(p+1));
        for (uint32_t i = _g_plen; i <= p; ++i)
        {
            mpz_init(_g_powers[i]);
            mpz_mul_ui(_g_powers[i],_g_powers[i-1],_g_base);
            _g_power_spmods[i] = (_g_power_spmods[i-1]*_g_base) % SPMOD;
        }
        _g_plen = p+1;
    }
    return _g_powers+p;
}

// returns base^p mod SPMOD
// must call get_power(p) first to ensure it is allocated
static inline spmod_t get_power_spmod(uint32_t p)
{
    return _g_power_spmods[p];
}

// grows the stack to ensure space for _g_stack[i]
static inline void ensure_stack_space(uint32_t i)
{
    if (i >= _g_slen)
    {
        _g_stack = realloc(_g_stack,sizeof(mpz_t)*(i+1));
        for (uint32_t j = _g_slen; j <= i; ++j)
            mpz_init(_g_stack[j]);
        _g_spmods = realloc(_g_spmods,sizeof(spmod_t)*(i+1));
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

// for temporary variables in functions
mpz_t _g_tmp0, _g_tmp1, _g_tmpU, _g_tmpV, _g_tmpQ;

/*
Dynamic memory handling for global variables. Call at the start and end.
*/

void init_globals()
{
    _g_buffer = malloc(BUFFER_SIZE);
    _g_buffer_index = 0;
    mpz_init_set_ui(MPZ_0,0);
    mpz_init_set_ui(MPZ_1,1);
    mpz_init_set_ui(MPZ_2,2);
    mpz_init_set_ui(MPZ_3,3);
    mpz_inits(_g_tmp0,_g_tmp1,_g_tmpU,_g_tmpV,_g_tmpQ,NULL);
    _g_stack = malloc(sizeof(mpz_t));
    mpz_init(_g_stack[0]);
    _g_spmods = malloc(sizeof(spmod_t));
    _g_slen = 1;
    _g_powers = malloc(sizeof(mpz_t));
    mpz_init_set_ui(_g_powers[0],1);
    _g_power_spmods = malloc(sizeof(spmod_t));
    _g_power_spmods[0] = 1;
    _g_plen = 1;
}

void clear_globals()
{
    free(_g_buffer);
    mpz_clears(MPZ_0,MPZ_1,MPZ_2,MPZ_3,NULL);
    mpz_clears(_g_tmp0,_g_tmp1,_g_tmpU,_g_tmpV,_g_tmpQ,NULL);
    for (uint32_t i = 0; i < _g_slen; ++i)
        mpz_clear(_g_stack[i]);
    free(_g_stack);
    free(_g_spmods);
    for (uint32_t i = 0; i < _g_plen; ++i)
        mpz_clear(_g_powers[i]);
    free(_g_powers);
    free(_g_power_spmods);
}

/*
Primality testing
*/

// PRP(2) test
// should avoid calling with even inputs for efficiency
// must not call with n == 2
static inline bool is_prime_prp2(const mpz_t n)
{
    mpz_sub_ui(_g_tmp0,n,1); // n-1
    mpz_powm(_g_tmp0,MPZ_2,_g_tmp0,n); // 2^(n-1) mod n
    return mpz_cmp_ui(_g_tmp0,1) == 0; // == 1
}

// SPRP(2) test
// must only call with odd n > 2
static inline bool is_prime_sprp2(const mpz_t n)
{
    mpz_sub_ui(_g_tmp0,n,1); // _g_tmp0 == n-1
    uint32_t s = mpz_scan1(_g_tmp0,0);
    mpz_fdiv_q_2exp(_g_tmp1,_g_tmp0,s); // n-1 == _g_tmp1 * 2^s
    mpz_powm(_g_tmp1,MPZ_2,_g_tmp1,n);
    if (mpz_cmp_ui(_g_tmp1,1) == 0 || mpz_cmp(_g_tmp1,_g_tmp0) == 0) // 1,n-1
        return true;
    while (--s) // s-1 squarings
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

// Strong Lucas test
// to be used as part of the BPSW test after the SPRP(2) test
static inline bool is_prime_lucas(const mpz_t n)
{
    int64_t D, Q;
    // iterate D = 5,-7,9,-11,... until (D|n) = -1 is found
    for (D = 5;;) // start at 5
    {
        // TODO optimize Jacobi symbol with spmod
        if (mpz_si_kronecker(D,n) < 0)
            break;
        D += 2; // next number (negative)
        if (mpz_si_kronecker(-D,n) < 0)
        {
            D = -D;
            break;
        }
        D += 2; // D = next number (positive)
        // if (D|n) = -1 not found quickkly, see if n is a square
        if (D == 129 && mpz_perfect_square_p(n))
            return false;
    }
    Q = (1-D)/4; // P = 1
    // factor 2s from n+1
    mpz_add_ui(_g_tmp0,n,1);
    uint32_t s = mpz_scan1(_g_tmp0,0);
    mpz_fdiv_q_2exp(_g_tmp0,_g_tmp0,s); // n+1 = _g_tmp0 * 2^s
    // initialize state for computing the _g_tmp0 index of U and V
    mpz_set_ui(_g_tmpU,1); // U_1 = 1
    mpz_set_ui(_g_tmpV,1); // V_1 = 1
    mpz_set_si(_g_tmpQ,Q); // Q^1 = Q
    mpz_mod(_g_tmpQ,_g_tmpQ,n); // handle negative Q (may not be necessary)
    uint32_t bit = mpz_sizeinbase(_g_tmp0,2)-1;
    while (bit) // TODO optimize bit testing
    {
        mpz_mul(_g_tmpU,_g_tmpU,_g_tmpV); // U_{2n} = U_n * V_n
        mpz_mod(_g_tmpU,_g_tmpU,n);
        mpz_powm_ui(_g_tmpV,_g_tmpV,2,n); // V_{n}^2
        mpz_submul_ui(_g_tmpV,_g_tmpQ,2); // V_{2n} = V_n^2 - 2*Q^n
        mpz_powm_ui(_g_tmpQ,_g_tmpQ,2,n); // Q^{2n} = (Q^n)^2
        if (mpz_tstbit(_g_tmp0,--bit)) // add 1 to the index
        {
            // store V_{n+1} in _g_tmp1 (not reduced modulo n yet)
            mpz_mul_si(_g_tmp1,_g_tmpU,D); // D*U_n
            mpz_add(_g_tmp1,_g_tmp1,_g_tmpV); // D*U_n + V_n
            if (mpz_odd_p(_g_tmp1)) // add n to make it even
                mpz_add(_g_tmp1,_g_tmp1,n);
            mpz_fdiv_q_2exp(_g_tmp1,_g_tmp1,1); // (D*U_n + V_n) / 2
            // store U_{n+1} in _g_tmpU
            mpz_add(_g_tmpU,_g_tmpU,_g_tmpV); // U_n + V_n
            if (mpz_odd_p(_g_tmpU)) // add n to make it even
                mpz_add(_g_tmpU,_g_tmpU,n);
            mpz_fdiv_q_2exp(_g_tmpU,_g_tmpU,1); // (U_n + V_n) / 2
            mpz_mod(_g_tmpU,_g_tmpU,n);
            // reduce computed V_{n+1} and store in _g_tmpV
            mpz_mod(_g_tmpV,_g_tmp1,n);
            // next power of Q
            mpz_mul_si(_g_tmpQ,_g_tmpQ,Q);
            mpz_mod(_g_tmpQ,_g_tmpQ,n);
        }
    }
    // with n+1 factored to d*2^s with d odd, perform the strong Lucas test
    if (mpz_cmp_ui(_g_tmpU,0) == 0)
        return true;
    if (mpz_cmp_ui(_g_tmpV,0) == 0)
        return true;
    while (--s) // s-1 doublings of V index
    {
        mpz_powm_ui(_g_tmpV,_g_tmpV,2,n); // V_n^2
        mpz_mul_2exp(_g_tmp1,_g_tmpQ,1); // 2*Q^n
        mpz_sub(_g_tmpV,_g_tmpV,_g_tmp1); // V_{2n} = V_n^2 - 2*Q^n
        mpz_mod(_g_tmpV,_g_tmpV,n);
        if (mpz_cmp_ui(_g_tmpV,0) == 0)
            return true;
        mpz_powm_ui(_g_tmpQ,_g_tmpQ,2,n); // Q^{2n} = (Q^n)^2
    }
    return false;
}

// Trial division test with small primes
static inline bool is_prime_tdiv(spmod_t spmod)
{
    return (spmod % 2)
        && (spmod % 3)
        && (spmod % 5)
        && (spmod % 7)
        && (spmod % 11)
        && (spmod % 13)
        && (spmod % 17)
        && (spmod % 19)
        && (spmod % 23)
        && (spmod % 29)
        && (spmod % 31)
        && (spmod % 37)
        && (spmod % 41)
        && (spmod % 43);
}

/*
Truncatable primes recursion functions
These write the subtree bytes and the end byte
The caller is responsible for writing the value byte(s)
*/

// convenience for accessing variables
#define STACK_CURR (_g_stack[_g_depth])
#define STACK_PREV (_g_stack[_g_depth-1])
#define SPMOD_CURR (_g_spmods[_g_depth])
#define SPMOD_PREV (_g_spmods[_g_depth-1])
#define POWER_CURR (*get_power(_g_rlen+_g_depth-1))
#define POWER_SPMOD_CURR (get_power_spmod(_g_rlen+_g_depth-1))

// primality test function to use
static inline bool prime_test(const mpz_t n, spmod_t spmod)
{
    if (mpz_fits_ushort_p(n))
    {
        if (spmod < 64) // test bit in 2^2 + 2^3 + ... + 2^61
            return 2891462833508853932uL & (1uL << (spmod));
        else if (spmod < TDIV_LIMIT)
            return is_prime_tdiv(spmod);
    }
    return is_prime_tdiv(spmod) && is_prime_sprp2(n) && is_prime_lucas(n);
}

// macro for all instances of calling the prime test function
#define PRIME_TEST_CURR prime_test(STACK_CURR,SPMOD_CURR)
//#define PRIME_TEST_CURR mpz_probab_prime_p(STACK_CURR,0)

// right truncatable (A024770 for base 10)
void primes_r()
{
    ++_g_depth;
    if (_g_depth <= _g_maxdepth)
    {
        ensure_stack_space(_g_depth);
        // left shift to create a 0 digit on the right
        mpz_mul_ui(STACK_CURR,STACK_PREV,_g_base);
        SPMOD_CURR = (SPMOD_PREV*_g_base) % SPMOD;
        for (uint32_t d = 1; d < _g_base; ++d)
        {
            // increment right digit
            mpz_add_ui(STACK_CURR,STACK_CURR,1);
            ++SPMOD_CURR;
            if (SPMOD_CURR == SPMOD)
                SPMOD_CURR = 0;
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
        ensure_stack_space(_g_depth);
        mpz_set(STACK_CURR,STACK_PREV);
        SPMOD_CURR = SPMOD_PREV;
        for (uint32_t d = 1; d < _g_base; ++d)
        {
            // increment left digit
            mpz_add(STACK_CURR,STACK_CURR,POWER_CURR);
            SPMOD_CURR = (SPMOD_CURR+POWER_SPMOD_CURR) % SPMOD;
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
        ensure_stack_space(_g_depth);
        // append left
        mpz_set(STACK_CURR,STACK_PREV);
        SPMOD_CURR = SPMOD_PREV;
        for (uint32_t d = 1; d < _g_base; ++d)
        {
            mpz_add(STACK_CURR,STACK_CURR,POWER_CURR);
            SPMOD_CURR = (SPMOD_CURR+POWER_SPMOD_CURR) % SPMOD;
            if (PRIME_TEST_CURR)
            {
                write_byte(0); // subtree
                write_byte(d);
                primes_lor();
            }
        }
        // append right
        mpz_mul_ui(STACK_CURR,STACK_PREV,_g_base);
        SPMOD_CURR = (SPMOD_PREV*_g_base) % SPMOD;
        for (uint32_t d = 1; d < _g_base; ++d)
        {
            mpz_add_ui(STACK_CURR,STACK_CURR,1);
            ++SPMOD_CURR;
            if (SPMOD_CURR == SPMOD)
                SPMOD_CURR = 0;
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

#define LAR_POWER_INDEX (_g_rlen + (_g_depth << 1) - 1)

// left and right truncatable (A077390 for base 10)
void primes_lar()
{
    ++_g_depth;
    // digit length is _g_rlen+2*_g_depth
    if ((_g_depth<<1) <= _g_maxdepth)
    {
        ensure_stack_space(_g_depth);
        mpz_mul_ui(STACK_CURR,STACK_PREV,_g_base); // shift left
        SPMOD_CURR = (SPMOD_PREV*_g_base) % SPMOD;
        for (uint32_t dl = 1; dl < _g_base; ++dl)
        {
            // increment left digit
            mpz_add(STACK_CURR,STACK_CURR,*get_power(LAR_POWER_INDEX));
            SPMOD_CURR = (SPMOD_CURR+get_power_spmod(LAR_POWER_INDEX)) % SPMOD;
            // right digit loop
            for (uint32_t dr = 1; dr < _g_base; ++dr)
            {
                mpz_add_ui(STACK_CURR,STACK_CURR,1);
                ++SPMOD_CURR;
                if (SPMOD_CURR == SPMOD)
                    SPMOD_CURR = 0;
                if (PRIME_TEST_CURR)
                {
                    write_byte(dl); // subtree
                    write_byte(dr);
                    primes_lar();
                }
            }
            // backtrack right digit increment
            mpz_sub_ui(STACK_CURR,STACK_CURR,_g_base-1);
            if (SPMOD_CURR >= _g_base-1)
                SPMOD_CURR -= _g_base-1;
            else
                SPMOD_CURR = SPMOD - (_g_base-1-SPMOD_CURR);
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
    write_byte(0); // root value
    if (byte2) // 2nd byte for root value
        write_byte(0);
    mpz_set_ui(_g_stack[0],root);
    _g_spmods[0] = root % SPMOD;
    _g_depth = 0;
    _g_rlen = 0;
    _g_maxdepth = (_g_maxlength >= _g_rlen) ? _g_maxlength - _g_rlen : 0;
    while (root) // initial digit length
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
        _g_spmods[0] = root % SPMOD;
        _g_depth = 0;
        _g_rlen = 1;
        _g_maxdepth = (_g_maxlength >= 1) ? _g_maxlength - 1 : 0;
        if (PRIME_TEST_CURR)
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
            _g_spmods[0] = (rootl*_g_base+rootr) % SPMOD;
            _g_depth = 0;
            _g_rlen = 2;
            _g_maxdepth = (_g_maxlength >= 2) ? _g_maxlength - 2 : 0;
            if (PRIME_TEST_CURR)
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
        write_byte(0); // root value
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
        write_byte(0); // root value
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
        write_byte(0); // root value
        write_byte(0);
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
        write_byte(0); // root value
        write_byte(0);
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
    _g_maxdepth = -1;
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
