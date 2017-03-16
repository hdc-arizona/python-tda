#!/usr/bin/env python

import functools

##############################################################################

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

    # must decide vertex order
    def __add__(self, other):
        return SimplicialComplex(self.face_set.union(other.face_set))
    def __mul__(self, other):
        return SimplicialComplex(self.face_set.intersection(other.face_set))

class Cover(object):
    def __init__(self, cpx, vertex_map):
        # build subcomplexes
        subcomplex_map = {}
        self.subcomplexes = {}
        self.complex = cpx
        for face in cpx.face_set:
            dest = min(vertex_map[i] for i in face)
            subcomplex_map.setdefault(dest, []).append(face)
        for (complex_id, subcomplex_faces) in subcomplex_map.items():
            self.subcomplexes[complex_id] = SimplicialComplex(subcomplex_faces)
        self.vertex_map = vertex_map

    # product here means the coarsest refinement of both covers
    # it's only defined when both covers cover the same simplicial complex
    def __mul__(self, other):
        new_vertex_map = {}
        for v in self.vertex_map.keys():
            new_vertex_map[v] = (self.vertex_map[v], other.vertex_map[v])
        return Cover(self.complex, new_vertex_map)

class BlowupComplex(object):

    def __init__(self, cover):
        self.cover = cover
        self.subcomplexes = {}
        
        # base subcomplexes
        for (complex_id, subcomplex) in cover.subcomplexes.items():
            self.subcomplexes[frozenset([complex_id])] = subcomplex

        base_complex_list = list(cover.subcomplexes.items())
        
        prev_set = self.subcomplexes
        
        # intersections

        while len(prev_set) > 0:
            next_set = {}
            # complex_id is a frozenset
            for (complex_id, subcomplex) in prev_set.items():
                # other_complex_id is a potential element of the sets in complex_id
                for (other_complex_id, other_subcomplex) in base_complex_list:
                    intersection_id = complex_id.union(frozenset([other_complex_id]))
                    if other_complex_id in complex_id:
                        # this base complex is already in our intersection, not interesting
                        continue
                    if intersection_id in next_set:
                        # we've already computed this, continue
                        continue
                    intersection_complex = subcomplex * other_subcomplex
                    if len(intersection_complex.face_set) > 0:
                        print("New intersection: %s" % intersection_id)
                        next_set[intersection_id] = intersection_complex
            self.subcomplexes.update(next_set)
            prev_set = next_set

        
##############################################################################

def from_faces(lst):
    return SimplicialComplex([frozenset(v) for v in lst])

def report_persistence(complex):
    red = complex.reduce_boundary_matrix()
    print("Reduced matrix")
    print(red[0])
    print("\nBetti numbers")
    for (b, i) in sorted(red[1].items()):
        print("  %d: %d" % (b, i))
    print("\nPersistence diagram")
    for (col, row) in sorted(red[2].items()):
        print("  Birth: %s - Death: %s (persistence: %d)" %
              (sorted(list(ex1.faces[row])),
               sorted(list(ex1.faces[col])), col - row))
    print(" Infinitely persistent components:")
    for i, col in enumerate(red[0]):
        if len(col) == 0 and i not in red[3]:
            print("  Birth: %s" % ex1.faces[i])

# print(ex1.boundary_matrix())


