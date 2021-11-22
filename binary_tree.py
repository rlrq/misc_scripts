## adapted from a solution in: https://stackoverflow.com/question/2598437/how-to-implement-a-binary-tree

class BinaryNode:
    def __init__(self, val):
        self.l = None
        self.r = None
        self.v = val
        ## self.comp outputs one of 3 strings: lt, gt, eq
        self.comp = lambda a, b: "eq" if a == b else "lt" if a < b else "gt"
    
    # this needs to be set by the user if not using values that can be directly compared
    def set_compare(self, f):
        self.comp = f
    
    def compare(self, other):
        return self.comp(self.v, other)
    
    def __lt__(self, other): return self.compare(other) == "lt"
    def __le__(self, other): return self.compare(other) in {"lt", "eq"}
    def __eq__(self, other): return self.compare(other) == "eq"
    def __ge__(self, other): return self.compare(other) in {"gt", "eq"}
    def __gt__(self, other): return self.compare(other) == "gt"

class BinaryTree:
    def __init__(self, compare):
        self.root = None
        self.compare = None
    
    def get_root(self):
        return self.root
    
    def add(self, val):
        if self.root is None:
            self.root = Node(val)
        else:
            self._add(val, self.root)
    
    def _add(self, val, node):
        if val == node:
            return
        elif val < node:
            if node.l is not None:
                self._add(val, node.l)
            else:
                node.l = BinaryNode(val)
        else:
            if node.r is not None:
                self._add(val, node.r)
            else:
                node.r = BinaryNode(val)
    
    def find(self, val, node):
        if self.root is not None:
            return self._find(val, self.root)
        else:
            return None
    
    def _find(self, val, node):
        if val == node:
            return node
        elif (val < node.v and node.l is not None):
            self._find(val, node.l)
        elif (val > node.v and node.r is not None):
            self._find(val, node.r)
    
    def delete_tree(self):
        ## let garbage collector do its thing
        self.root = None
    
    def print_tree(self):
        if self.root is not None:
            self._printTree(self.root)
    
    def _print_tree(self, node):
        if node is not None:
            self._print_tree(node.l)
            print(str(node.v) + ' ')
            self._print_tree(node.r)


