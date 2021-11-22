import itertools

## returns a in b
def within(a, b):
    a = [int(x) for x in a]
    b = [int(x) for x in b]
    return min(a) >= min(b) and max(a) <= max(b)

## returns a in any range in ranges
def within_any(a, ranges):
    for r in ranges:
        if within(a, r):
            return True
    return False

## convert ranges of [start, end) into set of positions
## e.g. [(1, 3), (6, 10)] --> {1, 2, 6, 7, 8, 9}
def ranges_to_pos(ranges):
    return set(itertools.chain(*[set(range(x[0], x[1])) for x in ranges]))

## convert set of positions into list of ranges
## e.g. {1, 2, 6, 7, 8, 9} --> [(1, 3), (6, 10)]
def pos_to_ranges(pos):
    if not pos: return []
    pos = sorted(pos)
    output = []
    start_v = pos[0]
    last_v = pos[0] - 1
    for i in range(len(pos)):
        v = pos[i]
        if v - last_v != 1:
            output.append((start_v, last_v + 1))
            start_v = v
        if i == len(pos) - 1:
            output.append((start_v, v + 1))
        last_v = v
    return output

## takes iterable of sets of pos and outputs a single set of pos
def pos_union(pos):
    return set(itertools.chain(*list(pos)))

def ranges_subtract(r1, r2):
    pos_subtract = ranges_to_pos(r1) - ranges_to_pos(r2)
    return pos_to_ranges(pos_subtract)

## e.g. [[(1, 3), (6, 9)], [(2, 3), (6, 10)]] --> {1, 2, 6, 7, 8, 9} --> [(1, 3), (6, 10)]
def ranges_union(ranges):
    pos_set = ranges_to_pos(itertools.chain(*ranges))
    return pos_to_ranges(pos_set)
