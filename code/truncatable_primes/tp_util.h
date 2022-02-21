#pragma once

#include <gmp.h>
#include <stdint.h>

// BPSW primality test
#define PRIME_TEST(n) mpz_probab_prime_p(n,0)

// check prime type
bool is_r_truncprime(const mpz_t a, uint32_t b);
bool is_l_truncprime(const mpz_t a, uint32_t b);
bool is_lor_truncprime(const mpz_t a, uint32_t b);
bool is_lar_truncprime(const mpz_t a, uint32_t b);

// stack frame
typedef struct
{
    mpz_t n; // the number to test for primality
    uint32_t i; // the next append
    uint32_t c; // number of child nodes
    // r/l use 1..base-1
    // lor uses 1..base-1 for left and base+1..2*base-1 for right
    // lar uses 1..base*base-1 for 0 root, base..base*base-1 otherwise
    // at 0, write root byte(s) and setup new value on the stack
    // at last value, write end byte and exit frame
    char v[2]; // root byte(s)
}
TP_FRAME;

// generator state (for nonzero roots only)
typedef struct
{
    // constants
    uint32_t base; // number system base >= 2, <= 65535
    mpz_t root; // root of the recursion
    uint32_t rlen; // length of root (digits in base)
    uint32_t maxlen; // maximum length of numbers allowed
    uint32_t mode; // 0 = bytes only, 1 = pre order, 2 = post order
    // expandable constants
    mpz_t *pow; // powers of the base (always initialized with length >= 1)
    uint32_t plen; // length of powers
    // recursion state
    uint32_t depth; // index in the stack
    TP_FRAME *stack; // number stack (always initialized with length >= 1)
    uint32_t slen; // length of stack
}
TP_STATE;

// yielded value from generator
// len == 0 for no value set, otherwise len != 0
typedef struct
{
    uint32_t len; // number of digits, 0 if value not yielded
    const mpz_t *num; // prime number, read only
    uint32_t children; // number of children (post order only)
    uint32_t path; // value for the append to parent node (TP_FRAME.i)
}
TP_VALUE;

#define TP_BYTES_ONLY 0
#define TP_PRE_ORDER 1
#define TP_POST_ORDER 2

// init/clear a TP_STATE
void TP_init(TP_STATE *state, uint32_t base, const mpz_t root,
        uint32_t maxlen, const char *rootv, uint32_t mode);
void TP_clear(TP_STATE *state);

// extract next element from generator
// state is the generator state, must be initialized/cleared
// ret is the next byte sequence, the return value is its length (1 or 2)
// value is the next number and info, not set on every call
// post == true for post order, enables info for number of child nodes
// return value is 0 at end of generator
uint32_t TP_next_r(TP_STATE *state, char *ret, TP_VALUE *value);
uint32_t TP_next_l(TP_STATE *state, char *ret, TP_VALUE *value);
uint32_t TP_next_lor(TP_STATE *state, char *ret, TP_VALUE *value);
uint32_t TP_next_lar(TP_STATE *state, char *ret, TP_VALUE *value);
