#pragma once

#include <stdint.h>
#include "functions.h"

// Requires preprocessor define of BASE

// If conditional, use as META_IF<uint64_t,condition,result_true,result_false>
// Does not work well recursively because C++ metaprogramming is not lazy
// This is not used below, due to the limitations with recursion.
template <class T, bool condition, T result_true, T result_false>
struct META_IF // true
{
    static constexpr T result = result_true;
};
template <class T, T result_true, T result_false>
struct META_IF<T,false,result_true,result_false> // false
{
    static constexpr T result = result_false;
};

/*
    Template metaprogramming for GCD. Use META_GCD<a,b>::result
*/
template <uint64_t a, uint64_t b> // modulus method
struct META_GCD
{
    static_assert(a != 0 || b != 0);
    static constexpr uint64_t result = META_GCD<b,a%b>::result;
};
template <uint64_t a>
struct META_GCD<a,0> // base case
{
    static_assert(a != 0);
    static constexpr uint64_t result = a;
};

/*
    Divides all factors of p from n. Use META_DIV_ALL<n,p>::result
*/
template <uint64_t n, uint64_t p, bool div> // div is (n%p)==0
struct META_DIV_ALL_H
{
    static constexpr uint64_t result =
        META_DIV_ALL_H<n/p,p,((n/p)%p)==0>::result;
};
template <uint64_t n, uint64_t p>
struct META_DIV_ALL_H<n,p,false> // n%p != 0
{
    static constexpr uint64_t result = n;
};
template <uint64_t n, uint64_t p>
struct META_DIV_ALL // abstract away the boolean variable from the client
{
    static_assert(p > 1);
    static constexpr uint64_t result = META_DIV_ALL_H<n,p,(n%p)==0>::result;
};

/*
    Product of distinct prime factors of n. This may need recursion depth up to
    O(sqrt(n)). Use META_DPF_PROD<n>::result
*/
template <uint64_t n, uint64_t d, uint64_t prod, bool in_range, bool div>
struct META_DPF_PROD_H // divisible, in range
{
    static_assert(d > 1);
    static_assert(n > 1);
    static constexpr uint64_t result =
        META_DPF_PROD_H<
            META_DIV_ALL<n,d>::result,
            d+1,
            prod*d,
            (d+1)*(d+1) <= META_DIV_ALL<n,d>::result,
            (n % (d+1)) == 0>::result;
};
template <uint64_t n, uint64_t d, uint64_t prod, bool div>
struct META_DPF_PROD_H<n,d,prod,false,div> // out of range, end
{
    static constexpr uint64_t result = prod*n;
};
template <uint64_t n, uint64_t d, uint64_t prod>
struct META_DPF_PROD_H<n,d,prod,true,false> // not divisible, increment div
{
    static constexpr uint64_t result =
        META_DPF_PROD_H<n,d+1,prod,((d+1)*(d+1)<=n),(n%(d+1))==0>::result;
};
template <uint64_t n>
struct META_DPF_PROD // abstract away other parameters from the client
{
    static_assert(n > 0);
    static constexpr uint64_t result =
        META_DPF_PROD_H<n,2,1,(4<=n),(n%2)==0>::result;
};

// Template metaprogramming to unroll coprime to base loop
// Use with META_LOOP<base>::function(n)
template <uint64_t i, uint64_t base, bool b> // line inclusion for i
struct META_LINE
{
    static inline void function(uint64_t n)
    {
        if (fermat_pp(n+i,base,mod_mult42)) printf("%lu\n",n+i);
    }
};
template <uint64_t i, uint64_t base>
struct META_LINE<i,base,false> // do not inline line for value i
{
    static inline void function(uint64_t){}
};
template <uint64_t i, uint64_t iters, uint64_t base, bool end>
struct META_LOOP_H
{
    static inline void function(uint64_t n)
    {
        META_LINE<i, base, META_GCD<base,i>::result == 1>::function(n);
        META_LOOP_H<i+1,iters,base,i+1==iters>::function(n);
    }
};
template <uint64_t i, uint64_t iters, uint64_t base>
struct META_LOOP_H<i,iters,base,true> // recursion base case
{
    static inline void function(uint64_t n){}
};
template <uint64_t base>
struct META_LOOP // abstract away boolean parameter from the client
{
    static inline void function(uint64_t n)
    {
        META_LOOP_H<0,META_DPF_PROD<base>::result,base,false>::function(n);
    }
};

