import numpy as np
from scipy.optimize import linear_sum_assignment
from fontTools.pens.recordingPen import RecordingPen
from fontTools.pens.statisticsPen import StatisticsPen

def contourStructure(contour):
    # Use second byte of the operation name (curveTo, closePath, etc),
    # as that's unique.
    return ''.join(op[0][1] for op in contour)

def outlineStructure(outline):
    return ''.join(contourStructure(contour) for contour in outline)

def outlinePosition(outline):
    assert(outline[0][0][0] == "moveTo")
    return outline[0][0][1][0]

def positionFlatOutline(outline, initPos):
    newOutline = []
    for op in outline:
        vec = []
        for x_y in op[1]:
            vec.append((x_y[0] + initPos[0], x_y[1] + initPos[1]))
        newOutline.append((op[0], tuple(vec)))
    return newOutline


def outlineVector(outline, flat=False):
    if not outline:
        return []
    if flat:
        outline = [outline]
    assert(outline[0][0][0] == "moveTo")
    initPos = outlinePosition(outline)
    vec = []
    for contour in outline:
        for op in contour:
            for x_y in op[1]:
                vec.append(x_y[0] - initPos[0])
                vec.append(x_y[1] - initPos[1])
    return tuple(vec)

def reconstructRecordingPenValues(structure, vector):
    # We saved the second char of the operation name; with num args
    op_map = {
      "o": ("moveTo", 1),
      "i": ("lineTo", 1),
      "C": ("qCurveTo", 2),
      "u": ("curveTo", 3),
      "l": ("closePath", 0),
      "n": ("endPath", 0),
    }
    ret = []
    it = iter(vector)
    for mnem in structure:
        op, nArgs = op_map[mnem]
        args = []
        for n in range(nArgs):
            x = next(it)
            y = next(it)
            args.append((x,y))
        ret.append((op, tuple(args)))
    try:
        next(it)
        assert False
    except StopIteration:
        pass
    return ret

def contourVector(c):
    rPen = RecordingPen()
    rPen.value = c
    stats = StatisticsPen()
    rPen.replay(stats)
    size = abs(stats.area) ** 0.5 * 0.5
    return (
        int(size),
        int(stats.meanX),
        int(stats.meanY),
        int(stats.stddevX * 2),
        int(stats.stddevY * 2),
        int(stats.correlation * size),
    )

def matchingCost(G, matching):
    return sum(G[i][j] for i, j in enumerate(matching))

def matchOutline(shape, ref, partial=False):
    if not partial:
        assert len(shape) == len(ref)
    if not len(shape): return shape, 0, []

    # Shortcut: If structures match assume it's correct.
    # Although if order is wrong matching can recover the
    # correct order...
    #if outlineStructure(shape) == outlineStructure(ref):
    #    return shape, 0, range(len(shape))

    # Perform a weighted-matching of outlines between shape and ref.
    # If found a perfect-matching, that's our solution.
    G = []
    refVecs = [np.array(contourVector(c)) for c in ref]
    shapeVecs = [np.array(contourVector(c)) for c in shape]
    refStructs = [contourStructure(c) for c in ref]
    shapeStructs = [contourStructure(c) for c in shape]

    for refStruct,refVec in zip(refStructs, refVecs):
        row = []
        G.append(row)
        for shapeStruct,shapeVec in zip(shapeStructs, shapeVecs):
            if refStruct != shapeStruct:
                row.append(1e10)
                continue
            diff = refVec - shapeVec
            diff = np.dot(diff, diff)
            row.append(diff)
    rows, cols = linear_sum_assignment(G)
    assert (rows == list(range(len(rows)))).all()
    cost = matchingCost(G, cols)
    if cost >= 1e10:
        return None, 1e10, None

    assignment = cols.tolist()

    # We have a matching. Reorder contours and return
    reordered = []
    matched = set()
    for c in cols:
        matched.add(c)
        reordered.append(shape[c])
    # Append any contours not matched, sorted by their structure
    other = []
    for i in range(len(shape)):
        if i in matched: continue
        other.append((shapeStructs[i], shape[i], i))
    other = sorted(other)
    for c in other:
        reordered.append(c[1])
        assignment.append(c[2])

    return reordered, cost, assignment

def reorderAssignment(lst, assignment):
    new = [None] * len(lst)
    for i,j in enumerate(assignment):
        new[i] = lst[j]
    return new
