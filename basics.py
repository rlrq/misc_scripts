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

def merge_2_overlapping_ranges(r1, r2):
    return min(r1 + r2), max(r1 + r2)

def merge_overlapping_ranges(*rs, **kwargs):
    last_len = len(rs)
    new_rs = []
    rs = list(rs)
    while True:
        while len(rs) > 0:
            r1 = rs.pop()
            i = 0
            while i < len(rs):
                if has_overlap(r1, rs[i], **kwargs):
                    r2 = rs.pop(i)
                    r1 = merge_2_overlapping_ranges(r1, r2)
                else:
                    i += 1
            new_rs.append(r1)
        if len(new_rs) == last_len:
            break
        else:
            last_len = len(new_rs)
            rs = new_rs
            new_rs = []
    return new_rs


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

def invert_dict(d):
    invd = {}
    for k, v in d.items():
        invd[v] = invd.get(v, []) + [k]
    return invd

# def merge_overlapping(*iters):
#     result = []
#     for s in iters:
#         s = set(s)
#         for t in result:
#             if t & s:
#                 t.update(s)
#                 break
#         else: ## no break
#             result.append(s)
#     return result

## returns list of sets
def merge_overlapping(*iters):
    last_len = len(iters)
    iters = [set(e) for e in iters]
    new_iters = []
    while True:
        while len(iters) > 0:
            s1 = iters.pop()
            i = 0
            while i < len(iters):
                if s1 & iters[i]:
                    s2 = iters.pop(i)
                    s1.update(s2)
                else:
                    i += 1
            new_iters.append(s1)
        if len(new_iters) == last_len:
            break
        else:
            last_len = len(new_iters)
            iters = new_iters
            new_iters = []
    return new_iters
