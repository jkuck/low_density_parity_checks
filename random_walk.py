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
    assert(np.abs(vector_count - nCr(n, w)) < .001), (vector_count, nCr(n, w)) 
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

#FIX CASE WHERE m != bin_count!!!

def plot_pr_Ax_zero(n, m, k, f, max_w):
    '''
    Inputs:
    - n: int, columns in the matrix A
    - m: int, rows in the matrix A
    - k: int, block size of 1's on the diagonal of A
    - f: float, probability with which remaining 0's are set to 1
    - w_max: int, plot Pr[A(x1-x2) = 0] up to (x1-x2) with hamming weight w_max

    '''
#    assert(n/k == m) #make more general later

    all_w = [] #w for the current partitioning
    all_pr_Ax_zero = [] #pr_Ax_zero for the current partitioning
    all_vector_count = [] #vectors in the current partitioning
    for w in range(1, max_w+1):
        print 'w =', w
        #(all_vc, Ax_zero_probs) = get_Ax_zero_probs(n, w, k, f, verbose=False)
        (all_vc, Ax_zero_probs) = get_Ax_zero_probs_incompleteCol(n, m, w, k, f, verbose=True)
        for prob, vector_count in Ax_zero_probs.iteritems():
            all_w.append(w)
            all_pr_Ax_zero.append(prob)
            all_vector_count.append(vector_count)

    #baseline where all elements of A are sampled with probability f + k/n
    #as in http://cs.stanford.edu/~ermon/papers/SparseHashing-revised.pdf
    f_prime = k/n + f
    baseline_w = []
    baseline_pr_Ax_zero = [] 
    baseline_vector_count = [] #vectors in the current partitioning

    assert(f_prime > f and f_prime <= .5), (f_prime, k, n, f)
    for w in range(1, max_w+1):
        baseline_w.append(w)
        prob = (.5 + .5*(1-2*f_prime)**w)**m
        baseline_pr_Ax_zero.append(prob)
        vector_count = nCr(n, w)
        baseline_vector_count.append(vector_count)

#    plt.scatter(all_w,all_pr_Ax_zero,c=all_vector_count)
#    print len(all_w)
#    print len(baseline_pr_Ax_zero)
#    print len(baseline_vector_count)
#    plt.scatter(baseline_w,baseline_pr_Ax_zero,c=baseline_vector_count,marker='x')
#    plt.colorbar()
#
#    plt.show()

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

    #calculate cumulative number of vectors with probability A(x1-x2)=0 greater than or equal to current value
    #note baseline is already monotonically decreasing with w, no need to sort
    cumulative_vector_count = 0
    baseline_cumulative_vector_counts = []
    for vc in baseline_vector_count:
        cumulative_vector_count += vc
        baseline_cumulative_vector_counts.append(cumulative_vector_count)

    print "baseline total vector count =", baseline_cumulative_vector_counts[-1]
    print "new total vector count =", cumulative_vector_counts[-1]
    print "2^n =", 2**n
    print "difference =", cumulative_vector_counts[-1] - baseline_cumulative_vector_counts[-1]
    print "2^n - (new total vector count)=", 2**n - cumulative_vector_counts[-1]
#    plt.scatter(all_w,all_pr_Ax_zero,c=cumulative_vector_counts)
#    print len(all_w)
#    print len(baseline_pr_Ax_zero)
#    print len(baseline_vector_count)
#    plt.scatter(baseline_w,baseline_pr_Ax_zero,c=baseline_cumulative_vector_counts,marker='x')
#    plt.colorbar()
#    plt.show()
#
    #calculate our upper bound on the sum of probabilities over a given set size
    probability_sum = 0
    probability_sums = []
    probability_averages = []
    for idx, vc in enumerate(all_vector_count):
        probability_sum += vc*all_pr_Ax_zero[idx]
        probability_sums.append(probability_sum)
        probability_averages.append(probability_sum/cumulative_vector_counts[idx])

    probability_sum = 0
    baseline_probability_sums = []
    baseline_probability_averages = []
    for idx, vc in enumerate(baseline_vector_count):
        probability_sum += vc*baseline_pr_Ax_zero[idx]
        baseline_probability_sums.append(probability_sum)
        baseline_probability_averages.append(probability_sum/baseline_cumulative_vector_counts[idx])


#    plt.scatter(cumulative_vector_counts, probability_sums, c = 'b', marker='+')
#    plt.scatter(baseline_cumulative_vector_counts, baseline_probability_sums, c='r', marker='x')
    plt.scatter(cumulative_vector_counts, probability_averages, c = 'b', marker='+')
    plt.scatter(baseline_cumulative_vector_counts, baseline_probability_averages, c='r', marker='x')
    plt.show()

def check_set_size_with_diag(n, m):
    print 



plot_pr_Ax_zero(n=40, m=4, k=1, f=.01, max_w=40)
#plot_pr_Ax_zero(n=120, m=30, k=4, f=.05, max_w=40)

count_vectors(n=20, w=5, k=4, verbose=True)
sleep(3)

partitions = integer_partitions(10)

print partitions
print type(partitions)
number_of_partitions = 0
for p in partitions:
    print p
    number_of_partitions += 1
print number_of_partitions
