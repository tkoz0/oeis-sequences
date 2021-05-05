/*

Program to enumerate probable primes. The numbers are written to stdout in base
10, 1 per line. Once the end is reached, it prints a line with "done".

TODO This can be better optimized for inputs larger than 42 bits.
TODO Base is limited to 32 bits to simplify some of the code for now (since it
will work with all the mod_mult functions). Expanding this to support 63 bit
bases is not very important for now.
TODO computing gcd() for optimizing is worse for primes or numbers with high
totient(n). This is faster by writing code that specifically only tests numbers
coprime to base without the need to compute gcd().

*/

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "functions.h"

static inline void loop_all(uint64_t min, uint64_t max, uint64_t base,
    uint64_t (*mod_mult)(uint64_t,uint64_t,uint64_t),
    bool (*test)(uint64_t,uint64_t,uint64_t(*)(uint64_t,uint64_t,uint64_t)))
{
    while (min <= max)
    {
        if (gcd32((uint32_t)base,(uint32_t)(min%base)) == 1 &&
            test(min,base,mod_mult))
            printf("%lu\n",min);
        ++min;
    }
    printf("done\n");
}

static inline void loop_odd(uint64_t min, uint64_t max, uint64_t base,
    uint64_t (*mod_mult)(uint64_t,uint64_t,uint64_t),
    bool (*test)(uint64_t,uint64_t,uint64_t(*)(uint64_t,uint64_t,uint64_t)))
{
    if (!(min&1)) ++min; // start at odd
    while (min <= max)
    {
        if (gcd32((uint32_t)base,(uint32_t)(min%base)) == 1 &&
            test(min,base,mod_mult))
            printf("%lu\n",min);
        min += 2;
    }
    printf("done\n");
}

void check_inputs(uint64_t min, uint64_t max, uint64_t base)
{
//    fprintf(stderr,"min=%lu\nmax=%lu\nbase=%lu\n",min,max,base);
    if (min < 2)
    {
        fprintf(stderr,"min number must be >= 2\n");
        exit(0);
    }
    if (max >= POW2(63))
    {
        fprintf(stderr,"max number is > 63 bits\n");
        exit(0);
    }
    if (base < 2)
    {
        fprintf(stderr,"base number must be >= 2\n");
        exit(0);
    }
    if (base >= POW2(32))
    {
        fprintf(stderr,"base number is > 32 bits\n");
        exit(0);
    }
    if (min > max)
    {
        fprintf(stderr,"must have min <= max\n");
        exit(0);
    }
}

int main(int argc, char **argv)
{
    if (argc < 5)
    {
        fprintf(stderr,"./a.out <min> <max> <base> <fpp|epp|ejpp|sfpp>\n");
        return 0;
    }
    // parse inputs
    uint64_t min = strtoul(argv[1],NULL,10);
    uint64_t max = strtoul(argv[2],NULL,10);
    uint64_t base = strtoul(argv[3],NULL,10);
    check_inputs(min,max,base);
//    fprintf(stderr,"type=%s\n",argv[4]);
    if (!strcmp("fpp",argv[4]))
    {
        if (max < POW2(32))
            loop_all(min,max,base,mod_mult32,fermat_pp);
        else if (max < POW2(42))
            loop_all(min,max,base,mod_mult42,fermat_pp);
        else
            loop_all(min,max,base,mod_mult63,fermat_pp);
    }
    else if (!strcmp("epp",argv[4]))
    {
        if (max < POW2(32))
            loop_odd(min,max,base,mod_mult32,euler_pp);
        else if (max < POW2(42))
            loop_odd(min,max,base,mod_mult42,euler_pp);
        else
            loop_odd(min,max,base,mod_mult63,euler_pp);
    }
    else if (!strcmp("ejpp",argv[4]))
    {
        if (max < POW2(32))
            loop_odd(min,max,base,mod_mult32,euler_jacobi_pp);
        else if (max < POW2(42))
            loop_odd(min,max,base,mod_mult42,euler_jacobi_pp);
        else
            loop_odd(min,max,base,mod_mult63,euler_jacobi_pp);
    }
    else if (!strcmp("sfpp",argv[4]))
    {
        if (max < POW2(32))
            loop_odd(min,max,base,mod_mult32,strong_fermat_pp);
        else if (max < POW2(42))
            loop_odd(min,max,base,mod_mult42,strong_fermat_pp);
        else
            loop_odd(min,max,base,mod_mult63,strong_fermat_pp);
    }
    else
        fprintf(stderr,"invalid type\n");
    return 0;
}
