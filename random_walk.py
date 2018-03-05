from __future__ import division
# Calculate Pr[A(x1-x2) = 0] when the (m x n) matrix A is sampled by 
#   1. beginning with a block diagonal matrix with k 1's in every row, 
#      e.g. for k = 4, m = 3, n = 12
#       [1 1 1 1 0 0 0 0 0 0 0 0]
#       [0 0 0 0 1 1 1 1 0 0 0 0]
#       [0 0 0 0 0 0 0 0 1 1 1 1]
#   2. setting remaining 0's to 1 with probability f

# When all elements are set to 1 with probability f, this probability can be
# computed using a biased random walk analysis (https://cs.stanford.edu/~ermon/papers/SparseHashing-revised.pdf).
# Here, we must consider the columns in which x1 and x2 differ, particular the parity of each k length chunk in
# x1 - x2.

import operator as op
import math
import matplotlib.pyplot as plt
import numpy as np
from bigfloat import BigFloat
import bigfloat as bf

def integer_partitions(n):
    '''
    Find integer partitions of n (http://jeromekelleher.net/generating-integer-partitions.html)

    Inputs:
    - n: type int, the integer to partition

    Outputs:
    - a generator where each item in the generator is a list representing a particular partition of n

    '''
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]

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

def get_Ax_zero_probs(n, w, k, f, verbose=False):
    '''
    Count the number of vectors of length n with hamming weight w by enumerating each way
    to split hamming weight between n/k bins of size k.

    Answer should be (n choose w) indpendent of k, double checking we get this answer using integer partitions

    Inputs:
    - n: int, the length of our binary vectors
    - w: int, hamming weight (number of ones in each vector)
    - k: int, bin size (see top), n should be a multiple of k
    - f: float, the probability that remaining elements in A are set to 1 

    Outputs:
    - vector_count: int, the number of binary vectors of length n with hamming weight w
    - Ax_zero_probs: dictionary with
        - key: float, probability for a particular partitioning of w that A(x1-x2) = 0
        - value: int, the number of vectors with hamming weight w in this particular partitioning
    '''
    bin_count = n//k
    assert(bin_count*k == n), (bin_count*k, n, bin_count, k)
    vector_count = 0
    weight_partitions = integer_partitions(w)

    Ax_zero_probs = {}
    for partition in weight_partitions:
        #we can split w into a maximum of bin_count separate bins
        #we can only put a maximum of k 1's in each bin
        if len(partition) <= bin_count and max(partition) <= k: 
            if verbose:
                print partition
            unique_summands = set(partition) #the unique summands in this particular partition

            cur_partition_bin_level_count = math.factorial(bin_count)
            for summand in unique_summands:
                cur_partition_bin_level_count /= math.factorial(partition.count(summand))
            cur_partition_bin_level_count /= math.factorial(bin_count - len(partition)) #order of empty bins also doesn't matter
            cur_partition_count = cur_partition_bin_level_count
            if verbose:
                print cur_partition_count            
            for summand in partition:
                cur_partition_count *= nCr(k, summand)
    
            #Calculate the probability that A(x1-x2) = 0 when (x1-x2) has hamming weight w
            #and w is partitioned among n/k bins as specified by partition
            pr_Ax_zero = 1.0 # the probability that A(x1-x2) = 0
            check_row_count = 0
            if verbose:
                print 
                print Ax_zero_probs
                print '-'*80
                print partition
            for summand in unique_summands:
                if verbose:
                    print "summand= ", summand
                row_prob = get_row_prob(summand, w, f)
                pr_Ax_zero *= row_prob**partition.count(summand)
                if verbose:
                    print "1row_prob =", row_prob
                    print "partition.count(summand) =", partition.count(summand)
                    print "row_prob**partition.count(summand) =", row_prob**partition.count(summand)
                check_row_count += partition.count(summand)
            row_prob = get_row_prob(0, w, f)
            pr_Ax_zero *= row_prob**(bin_count - len(partition))
            if verbose:
                print "2row_prob =", row_prob
                print "(bin_count - len(partition)) =", (bin_count - len(partition))
                print "row_prob**(bin_count - len(partition)) =", row_prob**(bin_count - len(partition))
    
                print "pr_Ax_zero =", pr_Ax_zero
            check_row_count += bin_count - len(partition)
            assert(check_row_count == bin_count), (check_row_count, bin_count) #bin_count is equal to the number of rows in A

            assert(not pr_Ax_zero in Ax_zero_probs)
            Ax_zero_probs[pr_Ax_zero] = cur_partition_count

            vector_count += cur_partition_count
            if verbose:
                print cur_partition_count

    # the total number of vectors should be n choose w, make sure this is correct
#    assert(np.abs(vector_count - nCr(n, w)) < .001), (vector_count, nCr(n, w)) 
    return (vector_count, Ax_zero_probs)


def get_Ax_zero_probs_incompleteCol(n, m, w, k, f, verbose=False):
    '''
    Count the number of vectors of length n with hamming weight w by enumerating each way
    to split hamming weight between m bins of size k where m < n/k. 

    Answer should be (n choose w) indpendent of k, double checking we get this answer using integer partitions

    Inputs:
    - n: int, the length of our binary vectors
    - w: int, hamming weight (number of ones in each vector)
    - k: int, bin size (see top), n should be a multiple of k
    - f: float, the probability that remaining elements in A are set to 1 

    Outputs:
    - vector_count: int, the number of binary vectors of length n with hamming weight w
    - Ax_zero_probs: dictionary with
        - key: float, probability for a particular partitioning of w that A(x1-x2) = 0
        - value: int, the number of vectors with hamming weight w in this particular partitioning
    '''
    #bin_count = n//k
    #assert(bin_count*k == n), (bin_count*k, n, bin_count, k)
    bin_count = m
    assert(bin_count < n/k)
    vector_count = 0

    Ax_zero_probs = {}

    for binned_weight in range(w+1):
        unbinned_weight = w - binned_weight #weight that is beyond 1's
        if binned_weight <= bin_count*k and unbinned_weight <= n - bin_count*k:

            bin_partitions = integer_partitions(binned_weight)

            for partition in bin_partitions:
                max_bin_weight = max(partition)
                if partition == [0]:
                    partition = [] #want partition of 0 to be empty list, not [0]
                #we can split w into a maximum of bin_count separate bins
                #we can only put a maximum of k 1's in each bin
                if len(partition) <= bin_count and max_bin_weight <= k: 
                    if verbose:
                        print partition
                    unique_summands = set(partition) #the unique summands in this particular partition

                    cur_partition_bin_level_count = math.factorial(bin_count)*nCr(n - bin_count*k, unbinned_weight)
                    for summand in unique_summands:
                        cur_partition_bin_level_count /= math.factorial(partition.count(summand))
                    cur_partition_bin_level_count /= math.factorial(bin_count - len(partition)) #order of empty bins also doesn't matter
                    cur_partition_count = cur_partition_bin_level_count
                    if verbose:
                        print cur_partition_count            
                    for summand in partition:
                        cur_partition_count *= nCr(k, summand)
            
                    #Calculate the probability that A(x1-x2) = 0 when (x1-x2) has hamming weight w
                    #and w is partitioned among n/k bins as specified by partition
                    pr_Ax_zero = 1.0 # the probability that A(x1-x2) = 0
                    check_row_count = 0
                    print 
                    print Ax_zero_probs
                    print '-'*80
                    print partition
                    for summand in unique_summands:
                        print "summand= ", summand
                        row_prob = get_row_prob(summand, w, f)
                        pr_Ax_zero *= row_prob**partition.count(summand)
                        print "1row_prob =", row_prob
                        print "partition.count(summand) =", partition.count(summand)
                        print "row_prob**partition.count(summand) =", row_prob**partition.count(summand)
                        check_row_count += partition.count(summand)
                    row_prob = get_row_prob(0, w, f)
                    pr_Ax_zero *= row_prob**(bin_count - len(partition))
                    print "2row_prob =", row_prob
                    print "(bin_count - len(partition)) =", (bin_count - len(partition))
                    print "row_prob**(bin_count - len(partition)) =", row_prob**(bin_count - len(partition))

                    print "pr_Ax_zero =", pr_Ax_zero
                    check_row_count += bin_count - len(partition)
                    assert(check_row_count == bin_count), (check_row_count, bin_count) #bin_count is equal to the number of rows in A

                    assert(not pr_Ax_zero in Ax_zero_probs)
                    Ax_zero_probs[pr_Ax_zero] = cur_partition_count

                    vector_count += cur_partition_count
                    if verbose:
                        print "cur_partition_count=", cur_partition_count

    # the total number of vectors should be n choose w, make sure this is correct
    assert(np.abs(vector_count - nCr(n, w)) < .001), (vector_count, nCr(n, w), vector_count-nCr(n, w)) 
    return (vector_count, Ax_zero_probs)


def get_row_prob(summand, w, f):
    '''
    Rows are sampled by setting a block of k elements to 1 and setting the remainder to 1 with probability f

    For a vector (x1 - x2) with hamming weight w and with summand 1's in the block set to 1, we find the probability
    that the dot product of (x1-x2) and the row is zero
    '''
    remaining_one_count = w - summand 
    if summand%2 == 0: #random walk starts at 0, want probability of ending up at 0 after remaining_one_count steps 
        row_prob = .5 + .5*(1-2*f)**remaining_one_count
    else: #random walk starts at 1, want probability of ending up at 0 after remaining_one_count steps 
        assert(summand%2 == 1)
        row_prob = .5 - .5*(1-2*f)**remaining_one_count

    #print "row_prob=", row_prob, "summand=", summand, "w=", w, "f=", f
    return row_prob


def bigFloat_nCr(n, r):
    '''
    Outputs:
    - ret_val: bigFloat, n choose r
    '''
    ret_val = bf.factorial(n)/(bf.factorial(r)*bf.factorial(n-r))
    return ret_val

def permutation_collision_prob_bigfloat(n, m, w_in, w_out, precision=100, debug=False):
    '''
    Define the set of matrices A_set with dimensions (m x n) with m <= n, exactly one 1 in each row,
    and zero or one 1's in each column.  Consider sampling a matrix A uniformly at random from A_set.

    Inputs:
    - n: int, columns in A
    - m: int, rows in A
    - w_in: int, hamming weight of a vector x which we will multiply with A
    - w_out: int, we're interested in whether the vector Ax has hamming weight w_out

    Outputs:
    - prob: BigFloat, the probability that Ax has hamming weight w_out for any vector x with hamming weight w_in
    '''
    assert(m <= n)
    assert(w_in <= n)
    assert(w_out >= m-(n-w_in))
    assert(w_out <= w_in and w_out <= m)
    with bf.precision(precision):
        #the number of matrices in A_set such that Ax has hamming weight w_out for any x with hamming weight w_in
        prob = bigFloat_nCr(w_in, w_out)*bigFloat_nCr(n-w_in, m-w_out) 
        #divide by the number of matrices in A_set
        prob /= bigFloat_nCr(n, m)
        if debug:
            A_set_size = bigFloat_nCr(n, m)
            check_A_set_size = 0
            for enum_w_out in range(max(0, m-(n-w_in)), min(w_in, m)+1):
                check_A_set_size += bigFloat_nCr(w_in, enum_w_out)*bigFloat_nCr(n-w_in, m-enum_w_out)
#            assert(np.abs(check_A_set_size - A_set_size) < .00001), (check_A_set_size, A_set_size)
        return prob

def permutation_collision_prob(n, m, w_in, w_out, return_ln=False, debug=False):
    '''
    Define the set of matrices A_set with dimensions (m x n) with m <= n, exactly one 1 in each row,
    and zero or one 1's in each column.  Consider sampling a matrix A uniformly at random from A_set.

    Inputs:
    - n: int, columns in A
    - m: int, rows in A
    - w_in: int, hamming weight of a vector x which we will multiply with A
    - w_out: int, we're interested in whether the vector Ax has hamming weight w_out
    - return_log: bool, if true return ln(prob) instead of prob

    Outputs:
    - prob: BigFloat, the probability that Ax has hamming weight w_out for any vector x with hamming weight w_in
    '''
    assert(m <= n)
    assert(w_in <= n)
    assert(w_out >= m-(n-w_in))
    assert(w_out <= w_in and w_out <= m)
    if debug:
        A_set_size = nCr(n, m)
        check_A_set_size = 0
        for enum_w_out in range(max(0, m-(n-w_in)), min(w_in, m)+1):
            check_A_set_size += nCr(w_in, enum_w_out)*nCr(n-w_in, m-enum_w_out)
        assert(check_A_set_size == A_set_size)
    if return_ln:
        ln = math.log
        #the number of matrices in A_set such that Ax has hamming weight w_out for any x with hamming weight w_in
        ln_prob = ln(nCr(w_in, w_out)) + ln(nCr(n-w_in, m-w_out))
        #divide by the number of matrices in A_set
        ln_prob -= ln(nCr(n, m))
        return ln_prob

    else:    
        #the number of matrices in A_set such that Ax has hamming weight w_out for any x with hamming weight w_in
        prob = nCr(w_in, w_out)*nCr(n-w_in, m-w_out) 
        #divide by the number of matrices in A_set
        prob /= nCr(n, m)
        return prob

#FIX CASE WHERE m != bin_count!!!
def plot_pr_Ax_zero(n, m, k, f_baseline, f_k1, max_w, BF_PRECISION=400, RUN_BASELINE=False):
    '''
    Inputs:
    - n: int, columns in the matrix A
    - m: int, rows in the matrix A
    - k: int, block size of 1's on the diagonal of A
    - f_baseline: float, probability with which 0's are set to 1 when running baseline with iid entries
    - f_k1: float, probability with which remaining 0's are set to 1 when k=1
    - w_max: int, plot Pr[A(x1-x2) = 0] up to (x1-x2) with hamming weight w_max

    '''
#    assert(n/k == m) #make more general later
    
    RUN_K_BLOCK_DIAG = True
    if RUN_K_BLOCK_DIAG:
        all_w = [1] #w for the current partitioning
        all_pr_Ax_zero = [1] #pr_Ax_zero for the current partitioning
        all_vector_count = [1] #vectors in the current partitioning
        for w in range(1, max_w+1):
            print 'w =', w
            #(all_vc, Ax_zero_probs) = get_Ax_zero_probs(n, w, k, f_baseline, verbose=False)
            (all_vc, Ax_zero_probs) = get_Ax_zero_probs_incompleteCol(n, m, w, k, f_baseline, verbose=False)
            for prob, vector_count in Ax_zero_probs.iteritems():
                all_w.append(w)
                all_pr_Ax_zero.append(prob)
                all_vector_count.append(vector_count)
        ###sort lists in order of probability A(x1-x2)=0
        all_pr_Ax_zero, all_vector_count, all_w = zip(*sorted(zip(all_pr_Ax_zero, all_vector_count, all_w)))
        #reverse lists so they are in order of descending probability A(x1-x2)=0
        all_pr_Ax_zero = list(reversed(all_pr_Ax_zero))
        all_vector_count = list(reversed(all_vector_count))
        all_w = list(reversed(all_w))
        #calculate cumulative number of vectors with probability A(x1-x2)=0 greater than or equal to current value
        cumulative_vector_count = 0
        cumulative_vector_counts = []
        for vc in all_vector_count:
            cumulative_vector_count += vc
            cumulative_vector_counts.append(cumulative_vector_count)

        #calculate our upper bound on the sum of probabilities over a given set size
        probability_sum = 0
        probability_sums = []
        probability_averages = []
        for idx, vc in enumerate(all_vector_count):
            probability_sum += vc*all_pr_Ax_zero[idx]
            probability_sums.append(probability_sum)
            probability_averages.append(probability_sum/cumulative_vector_counts[idx])


    if RUN_BASELINE:
        #baseline where all elements of A are sampled with probability f_baseline + k/n
        #as in http://cs.stanford.edu/~ermon/papers/SparseHashing-revised.pdf
        f_prime = k/n + f_baseline
    #    f_prime = f_baseline
        baseline_w = []
        baseline_pr_Ax_zero = [] 
        baseline_vector_count = [] #vectors in the current partitioning
    
        assert(f_prime > f_baseline and f_prime <= .5), (f_prime, k, n, f_baseline)
        for w in range(1, max_w+1):
            baseline_w.append(w)
            prob = (.5 + .5*(1-2*f_prime)**w)**m
            baseline_pr_Ax_zero.append(prob)
            vector_count = nCr(n, w)
            baseline_vector_count.append(vector_count)

    #calculate cumulative number of vectors with probability A(x1-x2)=0 greater than or equal to current value
    #note baseline is already monotonically decreasing with w, no need to sort
        cumulative_vector_count = 0
        baseline_cumulative_vector_counts = []
        for vc in baseline_vector_count:
            cumulative_vector_count += vc
            baseline_cumulative_vector_counts.append(cumulative_vector_count)

        probability_sum = 0
        baseline_probability_sums = []
        baseline_probability_averages = []
        for idx, vc in enumerate(baseline_vector_count):
            probability_sum += vc*baseline_pr_Ax_zero[idx]
            baseline_probability_sums.append(probability_sum)
            baseline_probability_averages.append(probability_sum/baseline_cumulative_vector_counts[idx])


    ####### 1. Begin with a block diagional with k = 1, so each row has one 1 and m columns have one 1
    ####### 2. Then randomly permute columns so that each row still has one 1 and columns still have one 1,
    ####### but the columns are random.
    ####### 3. With probability f change 0's to 1's
    permute_pr_Ax_zero = [] 
    permute_vector_count = [] #vectors in the current partitioning

    #NUMERICAL_METHOD = 'BigFloat'
    #NUMERICAL_METHOD = 'Logs'
    NUMERICAL_METHOD = 'Original'
    MULT_VECTOR_COUNT_IMPLICT = False
    if NUMERICAL_METHOD == 'BigFloat':
        for w in range(1, max_w+1):
            print 'w=', w
            total_vec_count = nCr(n, w)
            permute_vector_count.append(total_vec_count)
            with bf.precision(BF_PRECISION):
                calc_one = 0.0
                cur_prob_Ax_zero = 0
                for collision_count in range(max(0, m-(n-w)), min(w, m)+1):
                    if MULT_VECTOR_COUNT_IMPLICT:
                        calc_one = 1.0
                        cur_prob_Ax_zero += bigFloat_nCr(m, collision_count) * bigFloat_nCr(n - m, w - collision_count)\
                                            * ((.5 + .5*(1-2*f_k1)**w)**(m-collision_count)) * ((.5 - .5*(1-2*f_k1)**(w-1))**collision_count)

                    else:                    
                        collision_prob = permutation_collision_prob_bigfloat(n=n, m=m, w_in=w, w_out=collision_count, precision=BF_PRECISION, debug=True)
                        assert(collision_prob > 0), collision_prob
                        calc_one += collision_prob
                        cur_prob_Ax_zero += collision_prob * ((.5 + .5*(1-2*f_k1)**w)**(m-collision_count)) * ((.5 - .5*(1-2*f_k1)**(w-1))**collision_count)
                assert(np.abs(calc_one - 1.0) < .0001)            
                #assert(np.abs(cur_prob_Ax_zero - .5**m) < .001), (cur_prob_Ax_zero, .5**m)
                permute_pr_Ax_zero.append(cur_prob_Ax_zero)
                if (np.abs(total_vec_count*(.5**m - cur_prob_Ax_zero))>10**150):
                    print "w =", w, '!'*40
                    print 'cur_prob_Ax_zero = ', cur_prob_Ax_zero     
    elif NUMERICAL_METHOD == 'Logs':
        ln = math.log
        exp = math.exp
        for w in range(1, max_w+1):
            total_vec_count = nCr(n, w)
            permute_vector_count.append(total_vec_count)
            calc_one = 0.0
            cur_prob_Ax_zero = 0
            for collision_count in range(max(0, m-(n-w)), min(w, m)+1):
                ln_collision_prob = permutation_collision_prob(n=n, m=m, w_in=w, w_out=collision_count, return_ln=True, debug=True)
                assert(ln_collision_prob <= 0), ln_collision_prob
                calc_one += exp(ln_collision_prob)
                cur_prob_Ax_zero += ln_collision_prob * ((.5 + .5*(1-2*f_k1)**w)**(m-collision_count)) * ((.5 - .5*(1-2*f_k1)**(w-1))**collision_count)
            permute_pr_Ax_zero.append(cur_prob_Ax_zero)
            assert(np.abs(calc_one - 1.0) < .0001)
    else:
        assert(NUMERICAL_METHOD == 'Original')
        for w in range(1, max_w+1):
            total_vec_count = nCr(n, w)
            permute_vector_count.append(total_vec_count)
            USE_CORRECTED = False
            if USE_CORRECTED:
                calc_one = 0.0
                cur_prob_Ax_zero = 0
                for collision_count in range(max(0, m-(n-w)), min(w, m)+1):
                    collision_prob = permutation_collision_prob(n=n, m=m, w_in=w, w_out=collision_count, debug=True)
                    assert(collision_prob > 0), collision_prob
                    calc_one += collision_prob
                    cur_prob_Ax_zero += collision_prob * ((.5 + .5*(1-2*f_k1)**w)**(m-collision_count)) * ((.5 - .5*(1-2*f_k1)**(w-1))**collision_count)
#                assert(np.abs(cur_prob_Ax_zero - .5**m) < .001), (cur_prob_Ax_zero, .5**m)
                assert(np.abs(calc_one - 1.0) < .0001)
                permute_pr_Ax_zero.append(cur_prob_Ax_zero)
                if (np.abs(total_vec_count*(.5**m - cur_prob_Ax_zero))>10**150):
                    print "w =", w, '!'*40
                    print 'cur_prob_Ax_zero = ', cur_prob_Ax_zero
            else: #incorrect, was previously using, ACTUALLY equivalent
                double_check_vec_count = 0
                prob = 0
#                for i in range(min(w, m) + 1): # the number of elements in the vector (x1-x2) that hit a deterministic 1 in the matrix A
                for i in range(max(0, w-(n-m)), min(w, m) + 1):
                    #print "m=", m
                    #print "i=", i
                    #print "w=", w
                    #print 'n=', n
                    #print min(w, m, n-m) 
                    #if w-i > n-m:
                    #    continue
                    cur_vec_count = nCr(m, i)*nCr(n-m, w-i)
                    double_check_vec_count += cur_vec_count
                    cur_prob = cur_vec_count/total_vec_count
                    prob += cur_prob * ((.5 + .5*(1-2*f_k1)**w)**(m-i)) * ((.5 - .5*(1-2*f_k1)**(w-1))**i)
                assert(total_vec_count == double_check_vec_count)
                permute_pr_Ax_zero.append(prob)

    #calculate cumulative number of vectors with probability A(x1-x2)=0 greater than or equal to current value
    #note baseline is already monotonically decreasing with w, no need to sort
    cumulative_vector_count = 0
    permute_cumulative_vector_counts = []
    for vc in permute_vector_count:
        cumulative_vector_count += vc
        permute_cumulative_vector_counts.append(cumulative_vector_count)

    probability_sum = 0
    permute_probability_sums = []
    permute_probability_averages = []
    for idx, vc in enumerate(permute_vector_count):
        if MULT_VECTOR_COUNT_IMPLICT:
            probability_sum += permute_pr_Ax_zero[idx]
        else:
            probability_sum += vc*permute_pr_Ax_zero[idx]
        permute_probability_sums.append(probability_sum)
        permute_probability_averages.append(probability_sum/permute_cumulative_vector_counts[idx])


    #best possible bound, where f=.5
    #as in http://cs.stanford.edu/~ermon/papers/SparseHashing-revised.pdf
    f_best = .5
    best_w = []
    best_pr_Ax_zero = [] 
    best_vector_count = [] #vectors in the current partitioning

    for w in range(1, max_w+1):
        best_w.append(w)
        prob = (.5 + .5*(1-2*f_best)**w)**m
        best_pr_Ax_zero.append(prob)
        vector_count = nCr(n, w)
        best_vector_count.append(vector_count)
    #calculate cumulative number of vectors with probability A(x1-x2)=0 greater than or equal to current value
    #note best is already monotonically decreasing with w, no need to sort
    cumulative_vector_count = 0
    best_cumulative_vector_counts = []
    for vc in best_vector_count:
        cumulative_vector_count += vc
        best_cumulative_vector_counts.append(cumulative_vector_count)

    probability_sum = 0
    best_probability_sums = []
    best_probability_averages = []
    for idx, vc in enumerate(best_vector_count):
        probability_sum += vc*best_pr_Ax_zero[idx]
        best_probability_sums.append(probability_sum)
        best_probability_averages.append(probability_sum/best_cumulative_vector_counts[idx])




#    plt.scatter(all_w,all_pr_Ax_zero,c=all_vector_count)
#    print len(all_w)
#    print len(baseline_pr_Ax_zero)
#    print len(baseline_vector_count)
#    plt.scatter(baseline_w,baseline_pr_Ax_zero,c=baseline_vector_count,marker='x')
#    plt.colorbar()
#
#    plt.show()


    #calculate cumulative number of vectors with probability A(x1-x2)=0 greater than or equal to current value
    #note baseline is already monotonically decreasing with w, no need to sort
    if RUN_BASELINE and RUN_K_BLOCK_DIAG:  
        print "baseline total vector count =", baseline_cumulative_vector_counts[-1]
        print "new total vector count =", cumulative_vector_counts[-1]
        print "2^n =", 2**n
        print "difference =", cumulative_vector_counts[-1] - baseline_cumulative_vector_counts[-1]
        print "2^n - (new total vector count)=", 2**n - cumulative_vector_counts[-1]
        plt.scatter(all_w,all_pr_Ax_zero,c=cumulative_vector_counts)
        print len(all_w)
        print len(baseline_pr_Ax_zero)
        print len(baseline_vector_count)
        plt.scatter(baseline_w,baseline_pr_Ax_zero,c=baseline_cumulative_vector_counts,marker='x')
        plt.colorbar()
        plt.show()
        exit(0)

        print 'debugging'
        print len(permute_cumulative_vector_counts)
        print len(baseline_cumulative_vector_counts)

    for idx, permute_sum in enumerate(permute_probability_sums):
        print 'idx:', idx, 'permute_sum - best_probability_sums[idx]=', permute_sum - best_probability_sums[idx]
#        assert(permute_sum >= best_probability_sums[idx]), (idx, permute_sum, best_probability_sums[idx], (permute_sum - best_probability_sums[idx]))



#    plt.scatter(cumulative_vector_counts, probability_sums, c = 'b', marker='+')
    if RUN_BASELINE:
        print "len(baseline_cumulative_vector_counts):", len(baseline_cumulative_vector_counts)
        plt.scatter(baseline_cumulative_vector_counts, baseline_probability_sums, c='r', marker='x', label='baseline, iid f')
    plt.scatter(permute_cumulative_vector_counts, permute_probability_sums, c='g', marker='o', label='permuted k=1')
    plt.scatter(best_cumulative_vector_counts, best_probability_sums, c='b', marker='^', label='best, f=.5')
    plt.xlabel('set size')
    plt.ylabel('sum of p[Ax=0] for worst case set')
    plt.legend()
    plt.show()

    if RUN_BASELINE:
        print "len(baseline_cumulative_vector_counts):", len(baseline_cumulative_vector_counts)
        plt.scatter(baseline_cumulative_vector_counts, [baseline - best_probability_sums[idx] for idx, baseline in enumerate(baseline_probability_sums)], c='r', marker='x', label='baseline-best, iid f')
    plt.scatter(permute_cumulative_vector_counts, [permute - best_probability_sums[idx] for idx, permute in enumerate(permute_probability_sums)], c='g', marker='o', label='permuted-best k=1')
    #plt.scatter(best_cumulative_vector_counts, best_probability_sums, c='b', marker='^', label='best, f=.5')

#####    plt.scatter(cumulative_vector_counts, probability_averages, c = 'b', marker='+')
#####    plt.scatter(baseline_cumulative_vector_counts, baseline_probability_averages, c='r', marker='x')
#####    plt.scatter(permute_cumulative_vector_counts, permute_probability_averages, c='g', marker='o')
    plt.xlabel('set size')
    plt.ylabel('difference between methods of (sum of p[Ax=0] for worst case set)')
    
    plt.legend()
    plt.show()

def check_set_size_with_diag(n, m):
    print 

def quick_check_counts(n, m, w):
    sum1 = 0
    for i in range(w+1):
        sum1 += nCr(m, i)*nCr(n-m, w-i)

    print 'sum1 =', sum1
    print 'sum2 =', nCr(n, w)

#quick_check_counts(n=40, m=20, w=8)
#sleep(1)

#plot_pr_Ax_zero(n=2, m=2, k=1, f=.0, max_w=2)
if __name__=="__main__":
    plot_pr_Ax_zero(n=50, m=10, k=3, f_baseline=.01, f_k1=-1, max_w=50, BF_PRECISION=400, RUN_BASELINE=True)

    plot_pr_Ax_zero(n=300, m=80, k=1, f_baseline=.1, f_k1=.0966, max_w=300, RUN_BASELINE=True)
    exit(0)
    #check whether behavior with respect to n, m, f makes sense
    #n and max_w should be set equal
    plot_pr_Ax_zero(n=100, m=20, k=1, f_baseline=.01, f_k1=.01, max_w=100, RUN_BASELINE=True)
    plot_pr_Ax_zero(n=400, m=80, k=1, f_baseline=.01, f_k1=.01, max_w=400, RUN_BASELINE=True)
    plot_pr_Ax_zero(n=200, m=40, k=1, f_baseline=.01, f_k1=.01, max_w=200, RUN_BASELINE=True)

    plot_pr_Ax_zero(n=100, m=50, k=1, f_baseline=.02, f_k1=.02, max_w=100, RUN_BASELINE=True)        
    plot_pr_Ax_zero(n=576, m=5, k=1, f_baseline=.01, f_k1=.01, max_w=576, RUN_BASELINE=True)

###############    #plot_pr_Ax_zero(n=50, m=20, k=1, f_baseline=.1, f_k1=.5, max_w=50, RUN_BASELINE=False)
###############    plot_pr_Ax_zero(n=576, m=20, k=1, f_baseline=.4, f_k1=.4, max_w=576, RUN_BASELINE=True)
###############
###############    #plot_pr_Ax_zero(n=576, m=20, k=1, f_baseline=.2, f_k1=.2, max_w=576, RUN_BASELINE=True)
###############    plot_pr_Ax_zero(n=576, m=20, k=1, f_baseline=.1, f_k1=.05, max_w=576, RUN_BASELINE=True)
###############
###############    #why does this look messed up, can we do better than f=.5??
###############    plot_pr_Ax_zero(n=576, m=20, k=1, f_baseline=.1, f_k1=.1, max_w=576, RUN_BASELINE=True)
    
    ##count_vectors(n=20, w=5, k=4, verbose=True)
    ##sleep(3)
    ##
    ##partitions = integer_partitions(10)
    ##
    ##print partitions
    ##print type(partitions)
    ##number_of_partitions = 0
    ##for p in partitions:
    ##    print p
    ##    number_of_partitions += 1
    ##print number_of_partitions
    ##