/*
    Fermat Probable Prime finder, using template metaprogramming to compile
    efficient code for any base. The base is limited to [2,1023] to keep code
    size reasonable with the template metaprograming. Numbers are written to
    stdout, 1 per line, in base 10. This uses the 42 bit version for modular
    multiplication, capping to about 4.4 trillion. Once complete, a line
    containing "done" is written.
*/

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include "functions.h"
#include "templates.hpp"

// Probable prime base
#ifndef BASE
#error Must define base
#endif

// Range of supported bases
#if BASE < 2 || BASE > 1023
#error Base range allowed is [2,1023]
#endif

// Limit determined by choice of modular multiplication function
#define LIMIT (POW2(42)-1)

// Main loop in the mid range, assume midlo and midhi are divisible by BASE
// Each iteration handles integers coprime to BASE in (n,n+BASE)
static inline void loop_mid(uint64_t midlo, uint64_t midhi)
{
    uint64_t n;
    // Increment by distinct prime factor product as required for the loop
    // unrolling done by META_LOOP<>
    for (n = midlo; n != midhi; n += META_DPF_PROD<BASE>::result)
    {
        // Unroll a loop with length of the distinct prime factor product.
        // The numbers tested are coprime to the base.
        META_LOOP<BASE>::function(n);
    }
}

// Split main loop. Must be adjusted according to the base (by metaprogramming).
// Testing all of the few numbers at the endpoints barely affects performance.
static inline void loop(uint64_t min, uint64_t max)
{
    // Compute [midlo,midhi] range inside [min,max], endpoints divisible by base
    uint64_t midlo = min;
    if (midlo%BASE) midlo += BASE - (midlo%BASE);
    uint64_t midhi = max;
    midhi -= midhi%BASE;
    uint64_t n;
    // Below midlo, [min,midlo)
    for (n = min; n < midlo; ++n)
        if (fermat_pp(n,BASE,mod_mult42))
            printf("%lu\n",n);
    // Mid range, [midlo,midhi]
    loop_mid(midlo,midhi);
    // Above midhi, (midhi,max]
    for (n = midhi+1; n <= max; ++n)
        if (fermat_pp(n,BASE,mod_mult42))
            printf("%lu\n",n);
}

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr,"./a.out <min> <max>\n");
        return 0;
    }
    uint64_t min = strtoul(argv[1],NULL,10);
    uint64_t max = strtoul(argv[2],NULL,10);
    assert(2 <= min);
    assert(min <= max);
    assert(max <= LIMIT);
    loop(min,max);
    printf("done\n");
    return 0;
}
