#pragma once

#include <stdint.h>
#include <stdbool.h>

#define POW2(n) (1UL << (n))
#define BIT_MASK(n) (POW2(n) - 1)

/*
    Modular multiplication functions supporting varying bit sizes. Smaller
    modular multiplications can be done with fewer operations so it's best to
    use the smallest version sufficient.
*/

static inline uint64_t mod_mult32(uint64_t a, uint64_t b, uint64_t n)
{
    return (a*b)%n;
}

// split a and b into 21 bit parts to keep operations within 64 bits
static inline uint64_t mod_mult42(uint64_t a, uint64_t b, uint64_t n)
{
    uint64_t r = (a & BIT_MASK(21)) * b;
    r += (a >> 21) * ((b << 21) % n);
    return r % n;
}

// TODO add more versions for bit sizes between 42 and 63

// a full loop on individual bits
// a should be the smaller number to use fewer loop iterations
// this could probably be done faster with more than 64 bits
static inline uint64_t mod_mult63(uint64_t a, uint64_t b, uint64_t n)
{
    uint64_t r = 0;
    while (a)
    {
        if (a&1) r = (r+b)%n;
        // probably slower since this will result in about 2x divisions
        // r = (r + (b*(a&1))) % n;
        b = (b<<1)%n;
        a >>= 1;
    }
    return r;
}

/*
    GCD functions. Both inputs cannot be zero.
*/

static inline uint32_t gcd32(uint32_t a, uint32_t b)
{
    while (b)
    {
        uint32_t t = a%b;
        a = b;
        b = t;
    }
    return a;
}

static inline uint64_t gcd64(uint64_t a, uint64_t b)
{
    while (b)
    {
        uint64_t t = a%b;
        a = b;
        b = t;
    }
    return a;
}

/*
    Jacobi symbol calculator. Inputs are nonnegative and p must be odd > 1.
    Value returned will be in {-1,0,1}.
*/
static inline int8_t jacobi(uint64_t n, uint64_t p)
{
//    if (p == 1) return 1; // choosing not to handle p=1 case
    // computation state: compute (-1)^(sign-1) * (n|p)
    int8_t sign = 1; // sign is (-1)^(sign-1)
    uint64_t t;
    for (;;) // loop assumes p>1 is odd
    {
        n %= p; // n < p
        if (n == 0) return 0; // (0|p) = 0
        while (!(n&3)) n >>= 2; // divide out 4 since (4|p)=(2|p)^2=1
        if (!(n&1)) // even, must consider (2|p)
        {
            // (2|p) = 1 if p=+-1 (mod 8), -1 if p=+-3 (mod 8)
            // compute 1+(p%8), 2,4,6,8, 4 bit is set for -1 cases
            sign += ((1+(p&7))&4)>>2;
            n >>= 1;
        }
        // at this point, n will be odd
        if (n == 1) return ((sign&1)<<1)-1;
        // n and p are odd, swap using reciprocity
        sign += ((n&3) + (p&3) == 6);
        t = n;
        n = p;
        p = t;
    }
}

/*
    Probable primes tests. Currently only PRP ("weak" Fermat) and SPRP (strong
    Fermat). The mod_mult parameter should be chosen properly for the inputs.
    
    All Euler probable primes are Fermat probable primes.
    All Euler-Jacobi probable primes are Euler probable primes.
    All Strong probable primes are Euler-Jacobi probable primes.
    Strong -> Euler-Jacobi -> Euler -> Fermat
*/

/*
    Fermat Probable Prime Test
    Computes r = b^(n-1) % n
    If r != 1 then n is composite
    If r == 1 then n is probably prime
    The mod_mult parameter is a function pointer to support choosing different
    functions to efficiently perform modular multiplication depending on the bit
    sizes of n and b. The caller must guarantee that n and b are small enough.
    The caller must guarantee n >= 2. To improve performance, a helpful theorem
    is that if gcd(n,b) != 1 then b^(n-1) != 1 (mod n). For performance, avoid
    using the trivial bases 1 and 0.
*/
static inline bool fermat_pp(uint64_t n, uint64_t b,
    uint64_t (*mod_mult)(uint64_t,uint64_t,uint64_t))
{
    // exponent, mod base, result
    uint64_t e = n-1, mb = b%n, r;
    if (e&1) r = mb;
    else r = 1;
    while (e >>= 1) // r = b^(n-1)
    {
        mb = mod_mult(mb,mb,n);
        if (e&1) r = mod_mult(r,mb,n);
    }
    return r == 1;
}

/*
    Euler Probable Prime Test
    Caller must ensure n > 2 is odd. Similarly, for performance, avoid the
    trivial bases 0 and 1. Let n = 2q+1.
    b^(n-1) = b^(2q) = 1 (mod n) --> (b^q + 1)(b^q - 1) = 0 (mod n)
    This leads to the condition: b^q = 1 or -1 (mod n)
    Similarly, this can be improved by using the theorem stated before that
    gcd(n,b) != 1 -> b^(n-1) != 1 (mod n).
*/
static inline bool euler_pp(uint64_t n, uint64_t b,
    uint64_t (*mod_mult)(uint64_t,uint64_t,uint64_t))
{
    uint64_t e = n>>1, mb = b%n, r;
    if (e&1) r = mb;
    else r = 1;
    while (e >>= 1) // r = b^((n-1)/2)
    {
        mb = mod_mult(mb,mb,n);
        if (e&1) r = mod_mult(r,mb,n);
    }
    return r == 1 || r == n-1;
}

/*
    Euler-Jacobi Probable Prime Test
    Caller must ensure n > 2 is odd. Ignore trivial bases 0 and 1 for better
    performance. This is slightly stronger than the Euler probable prime test,
    returning b^((n-1)/2) = (b|n) (mod n) where (b|n) is the Jacobi symbol.
*/
static inline bool euler_jacobi_pp(uint64_t n, uint64_t b,
    uint64_t (*mod_mult)(uint64_t,uint64_t,uint64_t))
{
    uint64_t e = n>>1, mb = b%n, r;
    if (e&1) r = mb;
    else r = 1;
    while (e >>= 1) // r = b^((n-1)/2)
    {
        mb = mod_mult(mb,mb,n);
        if (e&1) r = mod_mult(r,mb,n);
    }
    int8_t j = jacobi(b,n);
    return (r == 1 && j == 1) || (r == n-1 && j == -1);
}

/*
    Strong Fermat Probable Prime Test
    Input n must be odd and n > 2, the caller must ensure this
    Does not make sense for even n because it is based on modular square roots
    Let n-1 = d * 2^s (d is odd, s >= 1 since n-1 is even)
    1. If b^d == 1 (mod n) then return true, otherwise
    2. If b^(d * 2^r) == -1 (mod n) for some 0 <= r < s, return true
    Similarly, the caller must ensure n and b are small enough for mod_mult.
    Trivial bases 1 and 0 should be skipped for performance reasons. The theorem
    gcd(n,b) != 1 implies b^(n-1) != 1 (mod n) is helpful because if n passes
    this test, it is also a (weak) Fermat probable prime.
*/
static inline bool strong_fermat_pp(uint64_t n, uint64_t b,
    uint64_t (*mod_mult)(uint64_t,uint64_t,uint64_t))
{
    uint32_t s = 0; // will count number of squarings needed
    uint64_t d = (n-1)>>1; // n-1 is even so do the first shift
    while (!(d&1)) ++s, d >>= 1; // compute s,d (actually s-1, counts squarings)
    // compute b^d (mod n)
    uint64_t mb = b%n;
    uint64_t r = mb; // start here instead of 1 since d is odd
    while (d >>= 1)
    {
        mb = mod_mult(mb,mb,n);
        if (d&1) r = mod_mult(r,mb,n);
    }
    // handles case 1 and case 2 with r = 0 (in description above)
    if (r == 1 || r == n-1) return true;
    while (s--) // finish case 2
        if ((r = mod_mult(r,r,n)) == n-1)
            return true;
    return false;
}
