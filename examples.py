from main import *

ex2 = from_faces([[1,2,3]])
ex3 = from_faces([[1,2,3],[2,3,4]])
ex4 = from_faces([[1,2,3], [1,2,4], [1,3,4], [2,3,4]])

def simple_examples():
    print("Complex 1:")
    report_persistence(ex1)
    print("\nComplex 2:")
    report_persistence(ex2)
    print("\nComplex 3:")
    report_persistence(ex3)
    print("\nComplex 4:")
    report_persistence(ex4)

##############################################################################


ex1 = from_faces([[1,2], [2,3],
                  [4,5], [5,6],
                  [7,8], [8,9],
                  [1,4], [4,7],
                  [2,5], [5,8],
                  [3,6], [6,9]])

c1 = Cover(ex1, {1: 1, 2: 2, 3: 2,
                 4: 1, 5: 2, 6: 2,
                 7: 1, 8: 2, 9: 2})

# print(c1.subcomplexes[1].faces)
# print(c1.subcomplexes[2].faces)

# report_persistence(c1.subcomplexes[1])
# print("\n\n\n")
# report_persistence(c1.subcomplexes[2])

# c2 = Cover(ex1, {1: 1, 2: 1, 3: 1,
#                  4: 2, 5: 2, 6: 2,
#                  7: 2, 8: 2, 9: 2})
# print("\n")
# print(c2.subcomplexes[1].faces)
# print(c1.subcomplexes[2].faces)
# c12 = c1 * c2
# print("\n")
# print("11: %s" % c12.subcomplexes[(1,1)].faces)
# print("12: %s" % c12.subcomplexes[(1,2)].faces)
# print("21: %s" % c12.subcomplexes[(2,1)].faces)
# print("22: %s" % c12.subcomplexes[(2,2)].faces)

bc1 = BlowupComplex(c1)

# for sc_id, sc in bc1.subcomplexes.items():
#     print("%s: %s" % (sc_id, sc.faces))

# print(bc1.reduce_boundary_matrix(frozenset([1])))
# print(bc1.reduce_boundary_matrix(frozenset([2])))
# print(bc1.reduce_boundary_matrix(frozenset([1,2])))

# print(bc1.subcomplexes[frozenset([1,2])].faces)

def compound_face(k):
    (k1, k2) = k
    return "([%s],[%s])" % (",".join(str(v) for v in k1),
                            ",".join(str(v) for v in k2))
def simple_face(k):
    return "[%s]" % (",".join(str(v) for v in k))

for face in bc1.faces:
    print("%s: %s" % (compound_face(face), bc1.face_id[face]))
print(bc1.reduce_boundary_matrix(frozenset([1])))

def report_sc(s):
    for face in bc1.subcomplexes[s].faces:
        print("%s: %s" % (simple_face(face), bc1.subcomplexes[s].face_id[face]))
    print("%s" % bc1.compute_mapping(s))

report_sc(frozenset([1]))
report_sc(frozenset([2]))
report_sc(frozenset([1,2]))

