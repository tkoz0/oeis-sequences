/* Utility functions */
num_length(b,n) =
{
    local(n_,l);
    n_ = n;
    l = 0;
    while(n_ > 0,
        l += 1;
        n_ = floor(n_/b);
    );
    return(l);
}

/* Truncatable prime verification functions */
isTP_R(b,n) =
{
    if (n < b,
        return(isprime(n)),
        return(isprime(n) && isTP_R(b,floor(n/b)))
    );
}
isTP_L(b,n) =
{
    if (n < b,
        return(isprime(n)),
        local(p,d);
        p = num_length(b,n)-1;
        d = floor(n/b^p);
        return(isprime(n) && isTP_L(b,n-d*b^p));
    );
}
isTP_LAR(b,n) =
{
    if (n < b*b,
        return(isprime(n)),
        local(p,d);
        p = num_length(b,n)-1;
        d = floor(n/b^p);
        return(isprime(n) && isTP_LAR(b,floor((n-d*b^p)/b)));
    );
}
isTP_LOR(b,n) =
{
    if (n < b,
        return(isprime(n)),
        local(p,d);
        p = num_length(b,n)-1;
        d = floor(n/b^p);
        return(isTP_LOR(b,floor(n/b)) || isTP_LOR(b,n-d*b^p));
    );
}

/*
Below functions get memory intensive and cannot be easily run with larger bases
*/

/* Perform all possible digit appends
b = base, d = length of n, n = prime */
appendTP_R(b,d,n) =
{
    /* n*b+r */
    local(ret,t);
    ret = List();
    for (r = 1, b-1,
        t = n*b+r;
        if (isprime(t), listinsert(ret,t,#ret+1));
    );
    return(ret);
}
appendTP_L(b,d,n) =
{
    /* l*b^d+n */
    local(ret,t);
    ret = List();
    for (l = 1, b-1,
        t = l*b^d+n;
        if (ispseudoprime(t), listinsert(ret,t,#ret+1));
    );
    return(ret);
}
appendTP_LAR(b,d,n) =
{
    /* l*b^(d+1)+n*b+r */
    local(ret,t);
    ret = List();
    for (l = 1, b-1,
        for (r = 1, b-1,
            t = l*b^(d+1)+n*b+r;
            if (isprime(t), listinsert(ret,t,#ret+1));
        );
    );
    return(ret);
}
appendTP_LOR(b,d,n) =
{
    /* l*b^d+n or n*b+r (duplicates included) */
    return(concat(appendTP_L(b,d,n),appendTP_R(b,d,n)));
}

/* BFS iteration, generates next layer */
tp_next_layer(b,d,nums,append_func) =
{
    local(ret);
    ret = List();
    for (i = 1, #nums,
        ret = concat(ret,append_func(b,d,nums[i]));
    );
    return(ret);
}

/* Creates list of the small primes */
list_1d_primes(b) = primes(primepi(b-1));
list_2d_primes(b) =
{
    local(tmp,start);
    tmp = list_1d_primes(b*b);
    start = #list_1d_primes(b)+1;
    return(tmp[start..#tmp]);
}

/* Generates all layers for a base */
gen_allTP_R(b) =
{
    local(ret,layer);
    ret = List();
    layer = list_1d_primes(b);
    while (#layer > 0,
        listinsert(ret,layer,#ret+1);
        layer = tp_next_layer(b,#ret,layer,appendTP_R);
    );
    return(ret);
}
gen_allTP_L(b) =
{
    local(ret,layer);
    ret = List();
    layer = list_1d_primes(b);
    while (#layer > 0,
        listinsert(ret,layer,#ret+1);
        layer = tp_next_layer(b,#ret,layer,appendTP_L);
    );
    return(ret);
}
gen_allTP_LAR(b) =
{
    local(ret,layer1,layer2);
    ret = List();
    layer1 = list_1d_primes(b);
    layer2 = list_2d_primes(b);
    while (#layer1 > 0 || #layer2 > 0,
        listinsert(ret,layer1,#ret+1);
        layer1 = tp_next_layer(b,#ret,layer1,appendTP_LAR);
        listinsert(ret,layer2,#ret+1);
        layer2 = tp_next_layer(b,#ret,layer2,appendTP_LAR);
    );
    return(ret);
}
gen_allTP_LOR(b,dups=1) =
{
    local(ret,layer);
    ret = List();
    layer = list_1d_primes(b);
    while (#layer > 0,
        listinsert(ret,layer,#ret+1);
        layer = tp_next_layer(b,#ret,layer,appendTP_LOR);
        if (!dups, layer = Set(layer));
    );
    return(ret);
}
