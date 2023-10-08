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


class Range(list):
    def __init__(self, start, end):
        super().__init__()
        self.start = min(start, end)
        self.end = max(start, end)
        self.extend([self.start, self.end])
    def __repr__(self):
        return f"Range({self.start}, {self.end})"
    ## convert to Ranges object
    def Ranges(self):
        return Ranges([self])
    def union(self, other):
        return self.Ranges().union(other.Ranges())
    def merge(self, other):
        if self.is_overlapping(other) or self.is_adjacent(other):
            return Range(min(self.start, other.start), max(self.end, other.end))
        else:
            raise Exception("Unable to merge disjoint ranges")
    def is_overlapping(self, other):
        return ((self.end > other.start and self.end <= other.end) or
                (self.start >= other.start and self.start < other.end))
    def is_adjacent(self, other):
        return (self.end == other.start or self.start == other.end)
    def offset(self, offset_val):
        return Range(self.start + offset_val, self.end + offset_val)
    def to_string(self, index = 0, end_exclusive = True):
        return f"{self.start + index}-{self.end + index - (0 if end_exclusive else 1)}"

## list of Range objects
class Ranges(list):
    ## ranges takes a list of Range objects. raw takes a list of lists and converts it list of Range objects.
    def __init__(self, ranges = [], raw = []):
        super().__init__()
        if raw:
            ranges = [Range(*r) for r in raw]
        self.extend(ranges)
        self.sort()
    def __repr__(self):
        return f"Ranges[{', '.join([str(r) for r in self])}]"
    def Ranges(self):
        return self
    def offset(self, offset_val):
        return Ranges([r.offset(offset_val) for r in self])
    def union(self, *others):
        all_ranges = self + list(itertools.chain(*others))
        all_ranges.sort()
        output = Ranges(all_ranges[:1])
        for r1 in all_ranges[1:]:
            tmp = Ranges()
            for i, r2 in enumerate(output):
                ## if r2 is before r1
                if r2.end < r1.start:
                    tmp.append(r2)
                ## else if r2 is after r1
                elif r2.start > r1.end:
                    tmp.append(r1)
                    tmp.extend(output[i:])
                    break
                ## else if overlapping
                else:
                    r1 = r1.merge(r2)
            tmp.append(r1)
            output = tmp
        return output
    def to_string(self, union = True, index = 0, end_exclusive = True):
        if union:
            ranges = self.union()
        else:
            ranges = Ranges(sorted(self))
        return ','.join([r.to_string(index = index, end_exclusive = end_exclusive) for r in ranges])
