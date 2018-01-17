import operator as op


def nCr(n, r):
    '''
    Compute n choose r
    https://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python
    '''
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom


def get_clause_count(XOR_var_count):
    parity_clause_count = 0
    for i in range(0, XOR_var_count+1, 2):
        parity_clause_count += nCr(XOR_var_count, i)
    print "parity_clause_count =", parity_clause_count
    print "2^XOR_var_count =", 2**XOR_var_count



get_clause_count(7)