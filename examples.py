from main import *

ex0 = from_faces([[1,2], [2,3]])
ex1 = from_faces([[1,2], [2,3],
                  [4,5], [5,6],
                  [7,8], [8,9],
                  [1,4], [4,7],
                  [2,5], [5,8],
                  [3,6], [6,9]])
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

# c1 = Cover(ex0, {1: 1, 2: 2, 3: 2})

c1 = Cover(ex1, {1: 1, 2: 2, 3: 2,
                 4: 1, 5: 2, 6: 2,
                 7: 1, 8: 2, 9: 2})

c2 = Cover(ex1, {1: 1, 2: 1, 3: 1,
                 4: 2, 5: 2, 6: 2,
                 7: 2, 8: 2, 9: 2})

c12 = c1 * c2

def compound_face(k):
    (k1, k2) = k
    return "([%s],[%s])" % (",".join(str(v) for v in k1),
                            ",".join(str(v) for v in k2))
def simple_face(k):
    return "[%s]" % (",".join(str(v) for v in k))

# bc1 = BlowupComplex(c1)

bc12 = BlowupComplex(c12)

print(bc12.reduce_boundary_matrix_partial([(1, 2), (2, 1), (2, 2)])[1])

