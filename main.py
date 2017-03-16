#!/usr/bin/env python

import functools

class Chain(object):

    def __init__(self, faces):
        if isinstance(faces, frozenset):
            self.d = faces
        elif isinstance(faces, Simplex):
            self.d = frozenset([faces.set])
        else:
            self.d = frozenset(faces)
    
    def boundary(self):
        result = Chain([])
        for f in self.d:
            result += Chain(Simplex(f).boundary())
        return result

    def __add__(self, other):
        return Chain(self.d.symmetric_difference(other.d))

    def __str__(self):
        return str(self.d)

# We basically never keep Simplex objects around, since they don't
# hash to their contents. Instead, we construct them, call their
# methods, and get the result.
class Simplex(object):

    def __init__(self, vals):
        if isinstance(vals, frozenset):
            self.set = vals
        else:
            self.set = frozenset(vals)

    def boundary(self):
        result = []
        for v in self.set:
            s = self.set - frozenset([v])
            if len(s) > 0:
                result.append(s)
        return result

##############################################################################

def set_lexicographical_cmp(s1, s2):
    if len(s1) < len(s2):
        return -1
    elif len(s1) > len(s2):
        return 1
    l1 = sorted(s1)
    l2 = sorted(s2)
    for el1, el2 in zip(l1, l2):
        if el1 < el2:
            return -1
        elif el1 > el2:
            return 1
    return 0

def reduce_column(left, right):
    result = []
    i = 0
    j = 0
    ll = len(left)
    lr = len(right)
    while i < ll and j < lr:
        if left[i] < right[j]:
            result.append(left[i])
            i += 1
        elif left[i] > right[j]:
            result.append(right[j])
            j += 1
        else:
            i += 1
            j += 1
    return result

##############################################################################

class SimplicialComplex(object):

    def __init__(self, faces, vertex_order=functools.cmp_to_key(set_lexicographical_cmp)):
        self.face_sorter = lambda k: (len(k), vertex_order(k))
        self.face_set = set(faces)
        # include all faces
        self.faces = []
        to_process = sorted(faces, key=lambda k: len(k))
        while len(to_process) > 0:
            old_to_process = to_process
            self.faces.extend(old_to_process)
            to_process = []
            for face in old_to_process:
                boundary = Simplex(face).boundary()
                for bf in boundary:
                    if bf not in self.face_set:
                        to_process.append(bf)
                        self.face_set.add(bf)
        self.face_set = frozenset(self.face_set)
        self.faces = sorted(self.face_set, key=self.face_sorter)
        self.face_id = {}
        for face in self.faces:
            self.face_id[face] = len(self.face_id)

    def boundary_matrix(self):
        cols = []
        for face in self.faces:
            b = [self.face_id[f] for f in Simplex(face).boundary()]
            cols.append(sorted(b))
        return cols

    def reduce_boundary_matrix(self):
        betti = {}
        bm = self.boundary_matrix()
        low     = {} # from col indices to low[i]
        low_inv = {} # from low[i] to col indices
        for i, col in enumerate(bm):
            if len(col) == 0:
                face = self.faces[i]
                dim = len(face)
                betti[dim-1] = betti.get(dim-1, 0) + 1
                continue
            mx = col[-1]
            while mx in low_inv:
                col_to_reduce = bm[low_inv[mx]]
                new_col = reduce_column(col_to_reduce, col)
                col = new_col
                bm[i] = col
                if len(col) == 0:
                    mx = None
                    break
                mx = col[-1]
            if mx is None:
                if i in low:
                    t = low[i]
                    del low_inv[t]
                    del low[i]
            else:
                low[i] = mx
                low_inv[mx] = i
            if len(col) == 0:
                face = self.faces[i]
                dim = len(face)
                betti[dim-1] = betti.get(dim-1, 0) + 1
            else:
                face = self.faces[i]
                dim = len(face)
                betti[dim-2] = betti.get(dim-2, 0) - 1
                
        return bm, betti, low, low_inv


##############################################################################

def from_faces(lst):
    return SimplicialComplex([frozenset(v) for v in lst])

def report_persistence(complex):
    red = complex.reduce_boundary_matrix()
    print(red[0])
    print(red[1])
    for (col, row) in sorted(red[2].items()):
        print("Birth: %s - Death: %s (persistence: %d)" % (ex1.faces[row], ex1.faces[col], col - row))
    print("Infinitely persistent components:")
    for i, col in enumerate(red[0]):
        if len(col) == 0 and i not in red[3]:
            print("Birth: %s" % ex1.faces[i])

ex1 = from_faces([[1,2], [2,3],
                  [4,5], [5,6],
                  [7,8], [8,9],
                  [1,4], [4,7],
                  [2,5], [5,8],
                  [3,6], [6,9]])
ex2 = from_faces([[1,2,3]])
ex3 = from_faces([[1,2,3],[2,3,4]])
ex4 = from_faces([[1,2,3], [1,2,4], [1,3,4], [2,3,4]])

print("Complex 1:")
report_persistence(ex1)
print("\nComplex 2:")
report_persistence(ex2)
print("\nComplex 3:")
report_persistence(ex3)
print("\nComplex 4:")
report_persistence(ex4)



# print(ex1.boundary_matrix())

