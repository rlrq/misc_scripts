import time

def time_elapsed(f):
    def helper(*args, **kwargs):
        start_time = time.time()
        output = f(*args, **kwargs)
        print("---", f.__name__, "--- %s seconds ---" % (time.time() - start_time))
        return output
    return helper
