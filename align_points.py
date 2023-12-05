import click
import numpy as np
import csv
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.spatial import KDTree

def kabsch_umeyama(A, B):
    assert A.shape == B.shape
    n, m = A.shape

    EA = np.mean(A, axis=0)
    EB = np.mean(B, axis=0)
    VarA = np.mean(np.linalg.norm(A - EA, axis=1) ** 2)

    H = ((A - EA).T @ (B - EB)) / n
    U, D, VT = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U) * np.linalg.det(VT))
    S = np.diag([1] * (m - 1) + [d])

    R = U @ S @ VT
    c = VarA / np.trace(np.diag(D) @ S)
    t = EA - c * R @ EB

    return R, c, t

def load_datafile(source):
    out = defaultdict(lambda: (list(),list()))
    with open(source) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        next(reader) # Ignore the header
        for x,y,circle,dataset in reader:
            dataset = dataset.strip()
            out[dataset][0].append((float(x),float(y)))
            out[dataset][1].append(circle.strip() == "True")

    return {dataset : (np.array(coords),np.array(circled)) for dataset, (coords, circled) in out.items()}

def align_with_fiducials(dataset, points):

    """ This code assumes the fiducials are oriented roughly as so:
    (b)

    (a)                 (c)

    Returns a tuple of the indices of a, b, and c, and a transformed copy of the points centered on the midpoint of a->c,
    and rotated so a->c is horizontal.
    """
    coords, circled = points
    fiducials = np.nonzero(circled)[0]
    if len(fiducials) != 3:
        raise Exception(f"Invalid number of fiducial markers in image {dataset}")
    a,b,c = fiducials[np.argsort(list(np.sum(np.linalg.norm(coords[fiducials] - coords[fiducials[i]], axis = 1)) for i in range(3)))]
    
    midpoint = 0.5 * (coords[a] + coords[c])
    xaxis = coords[c] - coords[a]
    xaxis /= np.linalg.norm(xaxis)
    m = np.vstack([xaxis,[-1 * xaxis[1], xaxis[0]]]).T
    return (a,b,c), (coords - midpoint) @ m

def rough_align(datasets):
    """ Performs a rough alignment stage for the datasets - using the fiducial markers, roughly
    align each dataset to the same position, and then match the point indices up across all datasets.
    Finally, return an initial guess for the true point position by averaging each group of points.
    """
    permutations = {}
    idx = None
    reference_set = None
    n = None
    output = None
    for dataset,points in datasets.items():
        _, transformed = align_with_fiducials(dataset,points)
        if idx is None:
            n = len(transformed)
            output = transformed.copy()
            reference_set = dataset
            permutations[dataset] = np.arange(n)
            idx = KDTree(transformed)
        else:
            hit = np.ones(n)
            permutation = np.zeros(n, dtype = int)

            for i, p in enumerate(transformed):
                # Find the closest point
                _, j = idx.query(p,1)
                # Make sure it hasn't been taken
                if not hit[j]:
                    print(f"Outlier detected in either {reference_set} or {dataset}")
                    return None
                hit[j] = 0
                permutation[j] = i
                
            permutations[dataset] = permutation
            output += transformed[permutation]

    output /= len(permutations)
    return output, permutations

def fine_align(datasets, previous, iterations = 10):
    n = len(datasets)
    for i in range(iterations):
        updated = np.zeros_like(previous)
        for points,_ in datasets.values():
            r,c,t = kabsch_umeyama(previous, points)
            updated += t + points @ r.T
        updated /= n
        delta = np.max(np.abs(previous - updated))
        previous = updated
        if delta < 1e-9:
            break

    return previous
    

_,ax = plt.subplots()
ax.set_aspect('equal')
    
data = load_datafile("cal_dat.csv")
# Use the fiducial markers to roughly align the datasets
average, permutations = rough_align(data)
average -= average.mean(axis = 0)
# Apply the permutations and forget the original order
for k,(points,circled) in data.items():
    p = permutations[k]
    data[k] = (points[p] - points.mean(axis = 0), circled[p])

# Align all of the points with a few iterations of aligning to an average, and then averaging the alignment
aligned = fine_align(data, average)
plt.scatter(aligned[:,0],aligned[:,1])
plt.show()
