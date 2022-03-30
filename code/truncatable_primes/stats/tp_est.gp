/* Truncatable primes estimation functions
R = right truncatable
L = left truncatable
LAR = left-and-right truncatable
LOR = left-or-right truncatable
*/

/* Estimated probability a d digit number in base b is prime
requires: b >= 2, d >= 1 */
P_unsimplified(b,d) = 1/d-1/((b-1)*d*(d-1));
P_approximated(b,d) = 1/d;
P(b,d) = P_approximated(b,d)/log(b);

/* Ratios for R,L,LAR */
ratio_R(b,d) = P(b,d)*b;
ratio_L(b,d) = P(b,d)*b*(b-1)/eulerphi(b);
ratio_LAR(b,d) = P(b,d)*b*(b-1);
/* Ratio for LOR, with duplicates included */
ratio_LOR_dup(b,d) = P(b,d)*b*(1+(b-1)/eulerphi(b));

/* Functions for estimating (without duplicates) next digit length */
est_C_R(b,d,prevC) = prevC*ratio_R(b,d);
est_C_L(b,d,prevC) = prevC*ratio_L(b,d);
est_C_LAR(b,d,prevC) = prevC*ratio_LAR(b,d);
est_C_LOR(b,d,prevC) = prevC*P(b,d)*b*((b-1)/eulerphi(b)+1-prevC/((b-1)^(d-1)*eulerphi(b)));
/* LOR estimation with duplicates */
est_C_LOR_dup(b,d,prevC) = prevC*ratio_LOR_dup(b,d);

/* Summations for R and L */
S(d,x) = (d!/x^d)*(exp(1)^x-sum(n=0,d,x^n/n!));
S_R(b,d) = S(d,b/log(b));
S_L(b,d) = S(d,(b*(b-1))/(log(b)*eulerphi(b)));

/* Summation for LOR with duplicates */
S_LOR_dup(b,d) = S(d,(b/log(b))*(1+(b-1)/eulerphi(b)));

/* Summation tests unsimplified formula direct implementation */
S_R_test(b,d,n) = sum(i=1,n,(b/log(b))^i/prod(j=1,i,d+j));
S_L_test(b,d,n) = sum(i=1,n,((b*(b-1))/(log(b)*eulerphi(b)))^i/prod(j=1,i,d+j));

/* default termination bound for numerical iteration */
eps_d = 0.001;

/* Estimate total number of d' digit truncatable primes for d' > d
b = base, d = num digits, C = value of C(b,d), eps = termination bound
returns a list, index 1 for d+1 digit, index 2 for d+2 digit ...
values continue as long as they are above the termination bound
*/
TP_est_helper(b,d,C,est_func,eps) =
{
    local(result,d_,Cprev,Cnext);
    result = List([]);
    d_ = d+1;
    Cprev = C;
    while ((Cnext = est_func(b,d_,Cprev)) > eps,
        listinsert(result,Cnext,#result+1);
        Cprev = Cnext;
        d_++;
    );
    return(result);
}
TP_est_all_R(b,d,C,eps=eps_d) = TP_est_helper(b,d,C,est_C_R,eps);
TP_est_all_L(b,d,C,eps=eps_d) = TP_est_helper(b,d,C,est_C_L,eps);
TP_est_all_LOR(b,d,C,eps=eps_d) = TP_est_helper(b,d,c,est_C_LOR,eps);
TP_est_all_LOR_dup(b,d,C,eps=eps_d) = TP_est_helper(b,d,c,est_C_LOR_dup,eps);
/* Cm1,C for C(b,d-1) and C(b,d)
continues as long as either of d-1+2k,d+2k are above termination bound*/
TP_est_all_LAR(b,d,Cm1,C,eps=eps_d) =
{
    local(result,d_,Cm1prev,Cm1next,Cprev,Cnext);
    result = List([]);
    d_ = d+2;
    Cm1prev = Cm1;
    Cprev = C;
    while (Cm1next = est_C_LAR(b,d_-1,Cm1prev); Cnext = est_C_LAR(b,d_,Cprev);
            Cm1next > eps || Cnext > eps,
        listinsert(result,Cm1next,#result+1);
        listinsert(result,Cnext,#result+1);
        Cm1prev = Cm1next;
        Cprev = Cnext;
        d_ += 2;
    );
    return(result);
}
