/*
Converts the tree output to integers, 1 per line, in the recursion order.

-i base (--input_base base)
    base of the truncatable primes (2-255, default 10)
-o base (--output_base base)
    base for text output (2-62, default 10, see GMP documentation)
-p type (--prime_type type)
    type of truncatable primes (r, l, lor, lar)
-r root (--root root)
    root number used for this recursion tree (default 0)
*/

#include <assert.h>
#include <ctype.h>
#include <getopt.h>
#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/*
Global constants and variables
*/

// command line arguments
const char *OPTION_STRING = "i:o:p:r:";
const struct option OPTION_LONG[] =
{
    { "input_base",  required_argument, 0, 'i' },
    { "output_base", required_argument, 0, 'o' },
    { "prime_type",  required_argument, 0, 'p' },
    { "root",        required_argument, 0, 'r' }
};

// buffer
#define BUFFER_SIZE (1<<16)
unsigned char *_g_buffer;
uint32_t _g_buffer_index, _g_buffer_bytes;

// read a byte from the buffer, -1 when EOF reached
static inline int read_byte()
{
    if (_g_buffer_index == _g_buffer_bytes)
    {
        _g_buffer_bytes = read(0,_g_buffer,BUFFER_SIZE);
        if (_g_buffer_bytes == -1)
        {
            fprintf(stderr,"unable to read input\n");
            exit(1);
        }
        _g_buffer_index = 0;
        if (_g_buffer_bytes == 0) // EOF
            return -1;
    }
    return _g_buffer[_g_buffer_index++];
}

// read a byte and error if reaching EOF
static inline int read_byte_strict()
{
    int b = read_byte();
    if (b == -1)
    {
        fprintf(stderr,"expected another byte but reoched EOF\n");
        exit(1);
    }
    return b;
}

// input arguments for bases and root
uint32_t _g_ibase = 10;
uint32_t _g_obase = 10;
uint64_t _g_root = 0;
char *_g_prime_type;

// for recursion
uint32_t _g_depth;
uint32_t _g_rlen;
mpz_t *_g_stack;
uint32_t _g_slen;
mpz_t *_g_powers;
uint32_t _g_plen;

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
            mpz_mul_ui(_g_powers[i],_g_powers[i-1],_g_ibase);
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

// macros for recursion usage
#define STACK_CURR (_g_stack[_g_depth])
#define STACK_PREV (_g_stack[_g_depth-1])
#define POWER_CURR (*get_power(_g_rlen+_g_depth-1))
#define CHECK_STACK ensure_stack_space(_g_depth)

// check bound of read byte, must be between (not inclusive)
static inline void check_byte_between(int b, int lo, int hi)
{
    if (b <= lo || b >= hi)
    {
        fprintf(stderr,"read byte out of bounds\n");
        exit(1);
    }
}

// write string of mpz_t to stdout
static inline void write_number(const mpz_t n)
{
    mpz_out_str(stdout,_g_obase,n);
    printf("\n");
}

/*
Truncatable prime tree reading functions.
Each reads the subtrees and the end byte.
*/

// right truncatable (A024770 for base 10)
void primes_r()
{
    int d = read_byte_strict(), dprev = 0;
    if (d == 255) return; // end
    ++_g_depth;
    CHECK_STACK;
    mpz_mul_ui(STACK_CURR,STACK_PREV,_g_ibase);
    do
    {
        check_byte_between(d,dprev,_g_ibase);
        mpz_add_ui(STACK_CURR,STACK_CURR,d-dprev);
        write_number(STACK_CURR);
        primes_r();
        dprev = d;
        d = read_byte_strict();
    }
    while (d != 255);
    --_g_depth;
}

// left truncatable (A024785 for base 10)
void primes_l()
{
    int d = read_byte_strict(), dprev = 0;
    if (d == 255) return; // end
    ++_g_depth;
    CHECK_STACK;
    mpz_set(STACK_CURR,STACK_PREV);
    do
    {
        check_byte_between(d,dprev,_g_ibase);
        mpz_addmul_ui(STACK_CURR,POWER_CURR,d-dprev);
        write_number(STACK_CURR);
        primes_l();
        dprev = d;
        d = read_byte_strict();
    }
    while (d != 255);
    --_g_depth;
}

// left or right truncatable (A137812 for base 10)
void primes_lor()
{
    int s = read_byte_strict(), sprev = -1;
    if (s == 255) return; // end
    int d, dprev = 0;
    ++_g_depth;
    CHECK_STACK;
    do
    {
        check_byte_between(s,-1,2);
        d = read_byte_strict(); // s=0 left, s=1 right, d=digit
        if (s == 0) // left
        {
            if (s != sprev) // switch append side
            {
                mpz_set(STACK_CURR,STACK_PREV);
                dprev = 0;
            }
            check_byte_between(d,dprev,_g_ibase);
            mpz_addmul_ui(STACK_CURR,POWER_CURR,d-dprev);
        }
        else // right
        {
            if (s != sprev) // switch append side
            {
                mpz_mul_ui(STACK_CURR,STACK_PREV,_g_ibase);
                dprev = 0;
            }
            check_byte_between(d,dprev,_g_ibase);
            mpz_add_ui(STACK_CURR,STACK_CURR,d-dprev);
        }
        write_number(STACK_CURR);
        primes_lor();
        dprev = d;
        sprev = s;
        s = read_byte_strict();
    }
    while (s != 255);
    --_g_depth;
}

#define LAR_POWER_INDEX (_g_rlen + (_g_depth << 1) - 1)

// left and right truncatable (A077390 for base 10)
void primes_lar()
{
    int ld = read_byte_strict(), ldprev = 0;
    if (ld == 255) return; // end
    int rd, rdprev = 0;
    ++_g_depth;
    CHECK_STACK;
    mpz_mul_ui(STACK_CURR,STACK_PREV,_g_ibase);
    do
    {
        rd = read_byte_strict();
        check_byte_between(ld,0,_g_ibase);
        if (ld == ldprev) // left digit is same (can be 0)
        {
            check_byte_between(rd,rdprev,_g_ibase);
            mpz_add_ui(STACK_CURR,STACK_CURR,rd-rdprev);
        }
        else // left digit changes
        {
            mpz_sub_ui(STACK_CURR,STACK_CURR,rdprev);
            check_byte_between(ld,ldprev,_g_ibase);
            mpz_addmul_ui(STACK_CURR,*get_power(LAR_POWER_INDEX),ld-ldprev);
            rdprev = 0;
            check_byte_between(rd,0,_g_ibase);
            mpz_add_ui(STACK_CURR,STACK_CURR,rd-rdprev);
        }
        write_number(STACK_CURR);
        primes_lar();
        ldprev = ld;
        rdprev = rd;
        ld = read_byte_strict();
    }
    while (ld != 255);
    --_g_depth;
}

void primes_lar_init()
{
    if (_g_rlen != 0)
        primes_lar();
    else // expect 1 or 2 digit roots, setup main recursion
    {
        int ld = read_byte_strict(), rd;
        if (ld == 255) return; // end
        int root, rootprev = 0;
        do
        {
            rd = read_byte_strict();
            check_byte_between(ld,-1,_g_ibase);
            check_byte_between(rd,-1,_g_ibase);
            root = ld*_g_ibase + rd;
            check_byte_between(root,rootprev,_g_ibase*_g_ibase);
            mpz_set_ui(_g_stack[0],root);
            _g_rlen = 1 + (ld != 0);
            _g_depth = 0;
            write_number(_g_stack[0]);
            primes_lar();
            rootprev = root;
            ld = read_byte_strict();
        }
        while (ld != 255);
    }
}

int main(int argc, char **argv)
{
    // setup
    _g_buffer = malloc(BUFFER_SIZE);
    _g_buffer_index = 0;
    _g_buffer_bytes = 0;
    _g_stack = malloc(sizeof(mpz_t));
    mpz_init(_g_stack[0]);
    _g_slen = 1;
    _g_powers = malloc(sizeof(mpz_t));
    mpz_init_set_ui(_g_powers[0],1);
    _g_plen = 1;
    // parse args
    if (argc < 2)
    {
        fprintf(stderr,"tree_convert <-p prime_type> "
                        "[-i input_base] [-o output_base] [-r root]\n");
        return 0;
    }
    int o;
    while ((o = getopt_long(argc,argv,OPTION_STRING,OPTION_LONG,NULL)) != -1)
    {
        switch(o)
        {
        case 'i': // input base
            if (!is_number(optarg))
            {
                fprintf(stderr,"input base must be a number\n");
                return 0;
            }
            _g_ibase = atoi(optarg);
            break;
        case 'o': // output base
            if (!is_number(optarg))
            {
                fprintf(stderr,"output base must be a number\n");
                return 0;
            }
            _g_obase = atoi(optarg);
            break;
        case 'p': // prime type
            _g_prime_type = optarg;
            break;
        case 'r': // root
            if (!is_number(optarg))
            {
                fprintf(stderr,"root must be a number\n");
                return 0;
            }
            _g_root = atoll(optarg);
            break;
        default:
            fprintf(stderr,"error parsing arguments\n");
            fprintf(stderr,"tree_convert <-p prime_type> "
                        "[-i input_base] [-o output_base] [-r root]\n");
            return 0;
        }
    }
    if (_g_ibase < 2 || _g_ibase > 255)
    {
        fprintf(stderr,"input base %u out of valid range (2-255)\n",_g_ibase);
        return 0;
    }
    if (_g_obase < 2 || _g_obase > 62)
    {
        fprintf(stderr,"output base %u out of valid range (2-62)\n",_g_obase);
        return 0;
    }
    // initialize recursion
    mpz_set_ui(_g_stack[0],_g_root);
    uint64_t root = _g_root;
    _g_rlen = 0;
    while (root)
        ++_g_rlen, root /= _g_ibase;
    _g_depth = 0;
    int b;
    b = read_byte();
    if (b != 255) // first root byte
    {
        fprintf(stderr,"invalid root byte %i, expected 255\n",b);
        return 1;
    }
    if (strcmp(_g_prime_type,"lor") == 0 || strcmp(_g_prime_type,"lar") == 0)
    {
        b = read_byte();
        if (b != 255) // second root byte
        {
            fprintf(stderr,"invalid root byte %i, expected 255\n",b);
            return 1;
        }
    }
    if (strcmp(_g_prime_type,"r") == 0)
        primes_r();
    else if (strcmp(_g_prime_type,"l") == 0)
        primes_l();
    else if (strcmp(_g_prime_type,"lor") == 0)
        primes_lor();
    else if (strcmp(_g_prime_type,"lar") == 0)
        primes_lar_init();
    else
    {
        fprintf(stderr,"invalid prime type: %s\n",_g_prime_type);
        return 0;
    }
    b = read_byte();
    if (b != -1) // should be at EOF
    {
        fprintf(stderr,"extra bytes found at end\n");
        return 1;
    }
    // clear
    free(_g_buffer);
    for (uint32_t i = 0; i < _g_slen; ++i)
        mpz_clear(_g_stack[i]);
    free(_g_stack);
    for (uint32_t i = 0; i < _g_plen; ++i)
        mpz_clear(_g_powers[i]);
    free(_g_powers);
    return 0;
}
