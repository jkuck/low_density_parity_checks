import multiprocessing
import time

def add_up(x):
    sum_ = 0
    i = 0
    while i < (x+30000000):
        sum_ += i
        i+=1
    return x

def sleep(x):
    time.sleep(3)
    return "done sleeping %d" % x

#pool = multiprocessing.Pool(4)
#pool = multiprocessing.Pool(multiprocessing.cpu_count())

THREADS = 16
pool = multiprocessing.Pool(THREADS)


results = [pool.apply_async( add_up, [i] ) for i in range(100)]

t0 = time.time()

for result in results:
    (sum_) = result.get()
#    print "type(sum_) =", type(sum_)
    print "result =", sum_

t1 = time.time()

print "threads =", THREADS
print "total time =", t1-t0
print "time/thread =", (t1-t0)/THREADS
