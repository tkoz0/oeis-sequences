/*
Truncatable primes utility functions
*/

#include <assert.h>
#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "tp_util.h"

/*
Prime type checker functions
*/

bool is_r_truncprime(const mpz_t a, uint32_t b)
{
    assert(b > 1);
    if (mpz_cmp_ui(a,0) <= 0)
        return false;
    mpz_t n;
    bool ret = true;
    mpz_init_set(n,a);
    while (mpz_cmp_ui(n,0) > 0) // test primality and truncate right
    {
        if (!PRIME_TEST(n) || mpz_div_ui(n,n,b) == 0)
        {
            ret = false;
            break;
        }
        mpz_div_ui(n,n,b);
    }
    mpz_clear(n);
    return ret;
}

bool is_l_truncprime(const mpz_t a, uint32_t b)
{
    assert(b > 1);
    if (mpz_cmp_ui(a,0) <= 0)
        return false;
    mpz_t n,d,p;
    bool ret = true;
    mpz_init_set(n,a);
    uint32_t l = 0;
    while (mpz_cmp_ui(n,0) > 0) // count digits
        ++l, mpz_div_ui(n,n,b);
    mpz_set(n,a);
    mpz_init_set_ui(d,b);
    mpz_init(p);
    mpz_pow_ui(p,d,l-1); // value of most significant digit
    while (mpz_cmp_ui(n,0) > 0) // truncate digits from left
    {
        // test primality (n) and ensure left digit (d) is nonzero
        if (!PRIME_TEST(n) || (mpz_div(d,n,p), mpz_cmp_ui(d,0) == 0))
        {
            ret = false;
            break;
        }
        mpz_submul_ui(n,p,mpz_get_ui(d)); // truncate
        mpz_div_ui(p,p,b); // value of next most significant digit
    }
    mpz_clears(n,d,p,NULL);
    return ret;
}

// recursive
bool is_lor_truncprime(const mpz_t a, uint32_t b)
{
    assert(b > 1);
    if (mpz_cmp_ui(a,b) < 0) // single digit
        return PRIME_TEST(a);
    if (!PRIME_TEST(a))
        return false;
    mpz_t tr,tl,d,p; // right/left truncations, division, power
    mpz_inits(tr,tl,d,p,NULL);
    uint32_t l = 0;
    mpz_set(tr,a);
    while (mpz_cmp_ui(tr,0) > 0) // count digits
    {
        ++l;
        if (mpz_div_ui(tr,tr,b) == 0) // digits must be nonzero
            return false;
    }
    // at start of each iteration
    // - tr = a prime to start with
    // - p = value of most significant digit
    mpz_set(tr,a);
    mpz_set_ui(d,b);
    mpz_pow_ui(p,d,l-1);
    bool ret;
    for (;;) // loop to avoid recursion where only 1 recursive call is made
    {
        if (mpz_cmp_ui(tr,b) < 0) // reached single digit prime
        {
            ret = true;
            break;
        }
        mpz_div(d,tr,p); // most significant digit
        mpz_set(tl,tr);
        mpz_submul_ui(tl,p,mpz_get_ui(d)); // set left truncated result
        mpz_div_ui(tr,tr,b); // set right truncated result
        mpz_div_ui(p,p,b);
        bool trp = PRIME_TEST(tr);
        bool tlp = PRIME_TEST(tl);
        if (trp)
        {
            if (tlp)
            {
                ret = is_lor_truncprime(tr,b) || is_lor_truncprime(tl,b);
                break;
            }
            else // use right truncation for next iteration
                continue;
        }
        else
        {
            if (tlp) // use left truncation for next iteration
            {
                mpz_set(tr,tl);
                continue;
            }
            else // cannot find parent prime
            {
                ret = false;
                break;
            }
        }
    }
    mpz_clears(tr,tl,d,p,NULL);
    return ret;
}

bool is_lar_truncprime(const mpz_t a, uint32_t b)
{
    assert(b > 1);
    if (mpz_cmp_ui(a,0) <= 0)
        return false;
    mpz_t n,d,p;
    bool ret = true;
    mpz_init_set(n,a);
    uint32_t l = 0;
    while (mpz_cmp_ui(n,0) > 0) // count digits
        ++l, mpz_div_ui(n,n,b);
    mpz_set(n,a);
    mpz_init_set_ui(d,b);
    mpz_init(p);
    mpz_pow_ui(p,d,l-1); // value of most significant digit
    while (mpz_cmp_ui(n,0) > 0)
    {
        // test primality (n) and ensure left digit (d) is nonzero
        if (!PRIME_TEST(n) || (mpz_div(d,n,p), mpz_cmp_ui(d,0) == 0))
        {
            ret = false;
            break;
        }
        mpz_submul_ui(n,p,mpz_get_ui(d)); // truncate left
        mpz_div_ui(n,n,b); // truncate right
        mpz_div_ui(p,p,b*b); // shift most significant digit value by 2
    }
    mpz_clears(n,d,p,NULL);
    return ret;
}

/*
Truncatable primes generator
base = is number system base
root = number for start of the recursion, nonzero only
maxlength = longest number length allowed, use -1 for unlimited
rootv = 2 bytes root value, only 1 used sometimes
mode = 0 for bytes only, 1 for pre order, 2 for post order
*/

void TP_init(TP_STATE *state, uint32_t base, const mpz_t root,
        uint32_t maxlen, const char *rootv, uint32_t mode)
{
    assert(2 <= base && base <= 255);
    state->base = base;
    assert(mpz_cmp_ui(root,0) > 0);
    mpz_init_set(state->root,root);
    state->maxlen = maxlen;
    state->mode = mode;
    // set rlen
    state->rlen = 0;
    mpz_t e;
    mpz_init_set(e,root);
    while (mpz_cmp_ui(e,0) > 0)
        mpz_div_ui(e,e,base), ++state->rlen;
    mpz_clear(e);
    // initialize powers array
    state->pow = malloc(sizeof(*state->pow));
    mpz_init_set_ui(state->pow[0],1);
    state->plen = 1;
    // initialize recursion state
    state->stack = malloc(2*sizeof(*state->stack));
    mpz_init_set(state->stack[0].n,root);
    mpz_init(state->stack[1].n);
    state->stack[0].i = 0;
    state->stack[1].i = 0;
    state->stack[0].v[0] = state->stack[0].v[1] = 0;
    state->stack[1].v[0] = rootv[0];
    state->stack[1].v[1] = rootv[1];
    state->stack[0].c = state->stack[1].c = 0;
    state->slen = 2;
    state->depth = 1; // start in next stack frame
}

void TP_clear(TP_STATE *state)
{
    uint32_t i;
    mpz_clear(state->root);
    for (i = 0; i < state->plen; ++i)
        mpz_clear(state->pow[i]);
    for (i = 0; i < state->slen; ++i)
        mpz_clear(state->stack[i].n);
    free(state->pow);
    free(state->stack);
}

// convenience macros
#define S_POW(i) (state->pow[i])
#define S_BASE (state->base)
#define S_CURR (state->stack[state->depth])
#define S_PREV (state->stack[state->depth-1])

void TP_resize_stack(TP_STATE *state, uint32_t len)
{
    assert(len > state->slen);
    state->stack = realloc(state->stack,len*sizeof(*state->stack));
    for (uint32_t i = state->slen; i < len; ++i)
    {
        mpz_init(state->stack[i].n);
        state->stack[i].i = 0;
        state->stack[i].v[0] = state->stack[i].v[1] = 0;
        state->stack[i].c = 0;
    }
    state->slen = len;
}

void TP_resize_powers(TP_STATE *state, uint32_t len)
{
    assert(len > state->plen);
    state->pow = realloc(state->pow,len*sizeof(*state->pow));
    for (uint32_t i = state->plen; i < len; ++i)
    {
        mpz_init(S_POW(i));
        mpz_mul_ui(S_POW(i),S_POW(i-1),S_BASE);
    }
    state->plen = len;
}

// sets value from generator conditionally
static inline void set_value(bool is_set, TP_VALUE *value, uint32_t len,
        const mpz_t *num, uint32_t children, uint32_t path)
{
    if (is_set)
    {
        value->len = len;
        value->num = num;
        value->children = children;
        value->path = path;
    }
    else
        value->len = 0;
}

uint32_t TP_next_r(TP_STATE *state, char *ret, TP_VALUE *value)
{
    if (state->depth == 0)
        return 0;
    for (;;) // perform steps until yielding next byte sequence
    {
        if (S_CURR.i == 0) // yield root
        {
            ++S_CURR.i;
            set_value(state->mode==TP_PRE_ORDER,value,
                state->rlen+state->depth-1,&S_PREV.n,-1,S_PREV.i-1);
            ret[0] = S_CURR.v[0];
            return 1;
        }
        else if (state->rlen + state->depth > state->maxlen)
            S_CURR.i = S_BASE; // skip deeper recursion
        else if (S_CURR.i < S_BASE) // digit append loop
        {
            if (S_CURR.i == 1) // setup
                mpz_mul_ui(S_CURR.n,S_PREV.n,S_BASE);
            mpz_add_ui(S_CURR.n,S_CURR.n,1);
            ++S_CURR.i;
            if (PRIME_TEST(S_CURR.n))
            {
                ++S_CURR.c;
                ++state->depth; // enter new stack frame
                if (state->depth >= state->slen)
                    TP_resize_stack(state,state->depth+1);
                S_CURR.v[0] = S_PREV.i-1; // digit
                S_CURR.i = S_CURR.c = 0;
            }
        }
        else // yield end
        {
            set_value(state->mode==TP_POST_ORDER,value,
                state->rlen+state->depth-1,&S_PREV.n,S_CURR.c,S_PREV.i-1);
            --state->depth;
            ret[0] = 255;
            return 1;
        }
    }
}

uint32_t TP_next_l(TP_STATE *state, char *ret, TP_VALUE *value)
{
    if (state->depth == 0)
        return 0;
    for (;;)
    {
        if (S_CURR.i == 0) // yield root
        {
            ++S_CURR.i;
            set_value(state->mode==TP_PRE_ORDER,value,
                state->rlen+state->depth-1,&S_PREV.n,-1,S_PREV.i-1);
            ret[0] = S_CURR.v[0];
            return 1;
        }
        else if (state->rlen + state->depth > state->maxlen)
            S_CURR.i = S_BASE; // skip deeper recursion
        else if (S_CURR.i < S_BASE) // digit append loop
        {
            if (S_CURR.i == 1) // setup
                mpz_set(S_CURR.n,S_PREV.n);
            if (state->rlen + state->depth > state->plen)
                TP_resize_powers(state,state->rlen+state->depth);
            mpz_add(S_CURR.n,S_CURR.n,state->pow[state->rlen+state->depth-1]);
            ++S_CURR.i;
            if (PRIME_TEST(S_CURR.n))
            {
                ++S_CURR.c;
                ++state->depth;
                if (state->depth >= state->slen)
                    TP_resize_stack(state,state->depth+1);
                S_CURR.v[0] = S_PREV.i-1; // digit
                S_CURR.i = S_CURR.c = 0;
            }
        }
        else // yield end
        {
            set_value(state->mode==TP_POST_ORDER,value,
                state->rlen+state->depth-1,&S_PREV.n,S_CURR.c,S_PREV.i-1);
            --state->depth;
            ret[0] = 255;
            return 1;
        }
    }
}

uint32_t TP_next_lor(TP_STATE *state, char *ret, TP_VALUE *value)
{
    if (state->depth == 0)
        return 0;
    for (;;)
    {
        if (S_CURR.i == 0) // yield root
        {
            ++S_CURR.i;
            set_value(state->mode==TP_PRE_ORDER,value,
                state->rlen+state->depth-1,&S_PREV.n,-1,S_PREV.i-1);
            ret[0] = S_CURR.v[0];
            ret[1] = S_CURR.v[1];
            return 2;
        }
        else if (state->rlen + state->depth > state->maxlen)
            S_CURR.i = S_BASE;
        else if (S_CURR.i < S_BASE) // append left
        {
            if (S_CURR.i == 1) // setup
                mpz_set(S_CURR.n,S_PREV.n);
            if (state->rlen + state->depth > state->plen)
                TP_resize_powers(state,state->rlen+state->depth);
            mpz_add(S_CURR.n,S_CURR.n,state->pow[state->rlen+state->depth-1]);
            ++S_CURR.i;
            if (PRIME_TEST(S_CURR.n))
            {
                ++S_CURR.c;
                ++state->depth;
                if (state->depth >= state->slen)
                    TP_resize_stack(state,state->depth+1);
                S_CURR.v[0] = 0;
                S_CURR.v[1] = S_PREV.i-1;
                S_CURR.i = S_CURR.c = 0;
            }
        }
        else if (S_CURR.i == S_BASE) // switch to append right
        {
            mpz_mul_ui(S_CURR.n,S_PREV.n,S_BASE);
            ++S_CURR.i;
        }
        else if (S_CURR.i < 2*S_BASE) // append right
        {
            mpz_add_ui(S_CURR.n,S_CURR.n,1);
            ++S_CURR.i;
            if (PRIME_TEST(S_CURR.n))
            {
                ++S_CURR.c;
                ++state->depth;
                if (state->depth >= state->slen)
                    TP_resize_stack(state,state->depth+1);
                S_CURR.v[0] = 1;
                S_CURR.v[1] = S_PREV.i-1-S_BASE;
                S_CURR.i = S_CURR.c = 0;
            }
        }
        else // yield end
        {
            set_value(state->mode==TP_POST_ORDER,value,
                state->rlen+state->depth-1,&S_PREV.n,S_CURR.c,S_PREV.i-1);
            --state->depth;
            ret[0] = 255;
            return 1;
        }
    }
}

uint32_t TP_next_lar(TP_STATE *state, char *ret, TP_VALUE *value)
{
    if (state->depth == 0)
        return 0;
    for (;;)
    {
        if (S_CURR.i == 0) // yield root
        {
            S_CURR.i = S_BASE; // skip so left append is nonzero
            set_value(state->mode==TP_PRE_ORDER,value,
                state->rlen+2*(state->depth-1),&S_PREV.n,-1,S_PREV.i-1);
            ret[0] = S_CURR.v[0];
            ret[1] = S_CURR.v[1];
            return 2;
        }
        else if (state->rlen + 2*state->depth > state->maxlen)
            S_CURR.i = S_BASE*S_BASE;
        else if (S_CURR.i < S_BASE*S_BASE)
        {
            if (S_CURR.i == S_BASE) // setup
                mpz_mul_ui(S_CURR.n,S_PREV.n,S_BASE);
            if (S_CURR.i % S_BASE == 0) // increment left digit
            {
                if (S_CURR.i != S_BASE) // reset right digit to 0
                    mpz_sub_ui(S_CURR.n,S_CURR.n,S_BASE-1);
                if (state->rlen + 2*state->depth > state->plen)
                    TP_resize_powers(state,state->rlen+2*state->depth);
                mpz_add(S_CURR.n,S_CURR.n,state->pow[state->rlen+2*state->depth-1]);
                ++S_CURR.i;
            }
            else // increment right digit, test primality
            {
                mpz_add_ui(S_CURR.n,S_CURR.n,1);
                ++S_CURR.i;
                if (PRIME_TEST(S_CURR.n))
                {
                    ++S_CURR.c;
                    ++state->depth;
                    if (state->depth >= state->slen)
                        TP_resize_stack(state,state->depth+1);
                    S_CURR.v[0] = (S_PREV.i-1)/S_BASE;
                    S_CURR.v[1] = (S_PREV.i-1)%S_BASE;
                    S_CURR.i = S_CURR.c = 0;
                }
            }
        }
        else // yield end
        {
            set_value(state->mode==TP_POST_ORDER,value,
                state->rlen+2*(state->depth-1),&S_PREV.n,S_CURR.c,S_PREV.i-1);
            --state->depth;
            ret[0] = 255;
            return 1;
        }
    }
}
