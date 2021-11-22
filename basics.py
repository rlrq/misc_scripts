def head(x, n=5):
    return tuple(x)[:n]

def tail(x, n=5):
    return tuple(x)[-n:]

def has_overlap(r1, r2, end_inclusive = False):
    r1 = sorted(r1)
    r2 = sorted(r2)
    if end_inclusive:
        return r2[0] <= r1[0] <= r2[1] or r1[0] <= r2[0] <= r1[1]
    else:
        return r2[0] <= r1[0] < r2[1] or r1[0] <= r2[0] < r1[1]

def any_overlap(r1, rs, **kwargs):
    helper_has_overlap = lambda r: has_overlap(r1, r)
    for r in rs:
        if has_overlap(r1, r):
            return True
    return False

## sorting key ('key' variable) should work on 'x' in 'for x in iterable'
def ties(iterable, sorting_function, tie_breaker = lambda x: x.items()[0], **kwargs):
    wanted_x = sorting_function(iterable, **kwargs)
    max_items = {k: v for k, v in iterable if key(k) == value} if type(iterable, dict) else \
                type(iterable)([x for x in iterable if sorting_key(x) == sorting_key(wanted_x)])
    return tie_breaker(max_items)

def print_iter(iterable):
    for entry in iterable:
        print(entry)
    return

def get_count_dict(iterable):
    output = {}
    for e in iterable:
        output[e] = output.get(e, 0) + 1
    return output

def extend_range(t, ext):
    """
    Extend a 2-element iterable of integers (decrease 1st element, increase 2nd), assuming ordered iterable.
    Returns tuple.
    """
    return (t[0] - ext, t[1] + ext)

def get_recursively(d, default, *keys):
    def helper(d, keys):
        key = keys[0]
        if key in d:
            if len(keys) == 1: return d[key]
            else: return helper(d[key], keys[1:])
        else:
            return default
    return helper(d, keys)
