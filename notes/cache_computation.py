try:
    import cPickle as pickle
except ImportError:
    import pickle

from var_computation import log_g, get_Ax_zero_log_probs

try:
    with open('.cache', 'rb') as cache_file:
        log_g.cache = pickle.load(cache_file)
except IOError:
    pass

tuples = [(352, 25), (939, 50), (1337, 27), (1413, 30)]
f = 0.05
for n, m in tuples:
    print(n, m)
    k = n // m
    t = [get_Ax_zero_log_probs(n, m, w, k, f) for w in range(1, n + 1)]
    if n // m > 3:
        k = 3
        t = [get_Ax_zero_log_probs(n, m, w, k, f) for w in range(1, n + 1)]
    with open('.cache', 'wb') as cache_file:
        pickle.dump(log_g.cache, cache_file)
