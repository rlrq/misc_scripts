## based on Ajami and Cohen (2019) Enumerating Minimal Weight Set Covers. Proceedings - International Conference on Data Engineering, 518-529

import copy

def make_examples():
    ## make example
    a, b, c, d, e = 'a', 'b', 'c', 'd', 'e'
    global s1, s2, s3, s4, S, U
    s1 = Set('s1', 2, [a, b, c])
    s2 = Set('s2', 3, [a, b, d])
    s3 = Set('s3', 2, [c, d])
    s4 = Set('s4', 5, [d, e])
    S = SetOfSets(s1, s2, s3, s4)
    U = S.elements
    return

def make_examples_2():
    a, b, c, d, e, f, g, h, i, j = range(10)
    global s1, s2, s3, s4, S, U
    s1 = Set('s1', 9, [a, b, c, d, e, f, g, h, i])
    s2 = Set('s2', 7, [a, b, c, d, e, f, g])
    s3 = Set('s3', 3, [h, i, j])
    s4 = Set('s4', 6, [e, f, g, h, i, j])
    S = SetOfSets(s1, s2, s3, s4)
    U = S.elements
    return

class e:
    ## a, b, c are all SetOfSets objects
    def __init__(self, a, b, c, tiebreaker = None):
        self.a = a
        self.b = b
        self.c = c
        self.tiebreaker = tiebreaker
    def __repr__(self):
        list_names = lambda setofsets:sorted(s.name for s in setofsets)
        return (f"{self.__class__.__name__}({list_names(self.a)},"
                f" {list_names(self.b)},"
                f" {list_names(self.c)})")
    @property
    def C(self): return self.a
    @property
    def Sfc(self): return self.b
    @property
    def So(self): return self.c
    @property
    def w(self): return self.C.w
    def copy(self):
        return self.__class__(self.a.copy(),
                              self.b.copy(),
                              self.c.copy(),
                              tiebreaker = self.tiebreaker)
    def _make_sorting_key(self):
        if self.tiebreaker:
            return lambda s: (s.w, self.tiebreaker(s))
        return lambda s: s.w
    ## dedicated functions so this class can be adapted to work when So is a SetOfSets of CollapsedNamedSets
    ## (or something like it that tracks equivalent sets and only removes one set of the CollapsedNamedSets
    ## but not the entire CollapsedNamedSets)
    def So_min(self):
        sorting_key = self._make_sorting_key()
        return min(self.So - self.C, key = sorting_key)
    def remove(self, *args, **kwargs):
        super().remove(*args, **kwargs)

class e1(e):
    def __init__(self, a, b, c):
        super().__init__(a, b, c)
    @property
    def Sf(self): return self.Sfc

class e2(e):
    def __init__(self, a, b, c):
        super().__init__(a, b, c)
    @property
    def Sc(self): return self.Sfc

class Q(list):
    def __init__(self):
        super().__init__()
    @property
    def top(self) -> e: return SetOfSets() if not self else self[0]
    def is_empty(self): return len(self) == 0
    def extract_min(self):
        return min(self, key = lambda s:s.w)
        # top_e = self.top
        # # sorting_key = self._make_sorting_key()
        # s = top_e.So_min()
        # new_Sf = type(top_e.So)(s)
        # new_So = top_e.So.copy()
        # new_So.remove(s)
        # new_C = top_e.C.copy()
        # new_C.add(s)
        # return type(top_e)(new_C, new_Sf, new_So)
    

## set class that tracks set name and set weight
class Set(set):
    def __init__(self, name, weight, elements):
        self.name = name
        self.weight = weight
        super().__init__(elements)
    def __repr__(self):
        return f"Set(name='{self.name}', {set(self)})"
        # return f"Set(name = {self.name}, weight = {self.weight}, length = {self.length})"
    ## hashable by name
    def __hash__(self):
        return hash(self.name)
    @property
    def w(self): return self.weight
    def copy(self):
        return Set(self.name,
                   self.weight,
                   list(self))

class SetOfSets(set):
    def __init__(self, *Sets):
        """
        Arguments:
            Sets (iter): iter of Set objects
        """
        # self._sets = {s.name: s for s in Sets}
        super().__init__(Sets)
        # super().__init__(set().union(*Sets))
    def __sub__(self, other):
        output = super().__sub__(other)
        return self.__class__(*output)
    @property
    def sets(self): return tuple(self)
    @property
    def elements(self): return set().union(*self)
    @property
    def set_names(self): return [s.name for s in self.sets]
    @property
    def w(self): return float("Inf") if not self else sum(s.w for s in self.sets)
    # def _update(self):
    #     """
    #     Regenerate's self's stored set values using elements of sets in self._sets
    #     """
    #     self -= self
    #     self.update(self.elements)
    # def remove(self, *Sets):
    #     for Set in Sets:
    #         if Set.name in self._sets and self._sets.get(Set.name) == Set:
    #             del self._sets[Set.name]
    #     ## regenerate self's stored set values using updated list of self.sets
    #     self._update()
    #     return
    # def add(self, *Sets):
    #     for Set in Sets:
    #         if Set.name in self._sets:
    #             print(f"Set name '{Set.name}' already in use")
    #         else:
    #             self._sets[Set.name] = Set
    #     ## regenerate self's stored set values using updated list of self.sets
    #     self._update()
    #     return
    def union(self, *args, **kwargs):
        output = super().union(*args, **kwargs)
        return self.__class__(*output)
    def intersection(self, *args, **kwargs):
        output = super().intersection(*args, **kwargs)
        return self.__class__(*output)
    def symmetric_difference(self, *args, **kwargs):
        output = super().symmetric_difference(*args, **kwargs)
        return self.__class__(*output)
    def remove(self, *Sets):
        self -= set(Sets)
        return
    def add(self, *Sets):
        self.update(Sets)
        return
    def copy(self):
        return self.__class__(*self.sets)

def wc_ratio(s1, s2):
    """
    Calculates ratio of weight of s1 to number of elements in s2 that can be covered by s1.
    (wc is shorthand for weight-cover, a name I came up with because the authors did not)
    """
    cover = s2.intersection(s1)
    if not cover: return float("Inf")
    return s1.w/len(cover)

def approx_min_SC(U, S):
    U = set(U)
    if set(U) - S.elements:
        return False
    C = SetOfSets()
    while U != set():
        s = min(S.sets, key = lambda s:wc_ratio(s, U))
        C.add(s)
        U -= s
    return C

def argmin(set_of_sets):
    return min(set_of_sets, key = lambda s:s.w)

# output = []
def enum_approx_order_SC(U, S, set_of_sets = SetOfSets, num_sc = 100):
    # global output
    output = []
    Q1, Q2 = Q(), Q()
    C = approx_min_SC(U, S)
    if C:
        Q1.append(e1(C, set_of_sets(), S))
    ## limit number of iterations
    for i in range(100):
        ## exit if queues are empty
        if (Q1.is_empty() and Q2.is_empty()): break
        ## if first entry of Q1 has lower (or equal) C weight to Q2's
        if Q1.top.w <= Q2.top.w:
            ## get entry with lowest C weight
            top_e = Q1.extract_min()
            C, Sf, So = top_e.C, top_e.Sf, top_e.So
            output.append(C)
            ## some algorithm things that I don't really understand
            if (So - C) != set():
                s = argmin(So - C)
                new_C = C.union({s})
                new_Sf = set_of_sets(s)
                new_So = set_of_sets(*(So - C - new_Sf))
                Q2.append(e2(new_C, new_Sf, new_So))
            for s in C - Sf:
                So = So - set_of_sets(s)
                U_prime = U - Sf.elements
                C_prime = approx_min_SC(U_prime, So)
                if C_prime is not False:
                    Q1.append(e1(C_prime.union(Sf), Sf, So))
                Sf = Sf.union({s})
            null = Q1.pop(0)
        ## else if first entry of Q2 has lower weight C
        else:
            ## get entry with lowest C weight
            top_e = Q2.extract_min()
            C, Sc, So = top_e.C, top_e.Sc, top_e.So
            output.append(C)
            if So != set():
                s_prime = argmin(So)
                new_C = C.union({s_prime})
                new_Sc = set_of_sets(s_prime)
                Q2.append(e2(new_C - Sc, new_Sc.copy(), So - new_Sc))
                Q2.append(e2(new_C, new_Sc.copy(), So - new_Sc))
            null = Q2.pop(0)
    return output

make_examples_2()
x = enum_approx_order_SC(U, S)
[(setofsets.w, list(s.name for s in setofsets)) for setofsets in x[:10]]
## now we just gotta find a way to tie break this
