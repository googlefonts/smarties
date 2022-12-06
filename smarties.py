
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import RecordingPen
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.svgPathPen import main as svgMain
from fontTools.varLib.interpolatable import PerContourPen
from collections import Counter, defaultdict
import numpy as np
import math
from pprint import pprint
import sys

fontfile = sys.argv[1]
count = None
demoS = None
if len(sys.argv) > 2:
    arg = sys.argv[2]
    if arg.startswith("U+"):
        demoS = int(arg[2:], 16)
    else:
        count = int(arg)

font = TTFont(fontfile)
upem = font['head'].unitsPerEm
hhea = font['hhea']
ascent, descent = hhea.ascent, hhea.descent
cmap = font['cmap'].getBestCmap()

LBase = 0x1100
VBase = 0x1161
TBase = 0x11A7
LCount = 19
VCount = 21
TCount = 28
SBase = 0xAC00
NCount = VCount * TCount
SCount = count if count else LCount * NCount

def isCombiningL(u):
    return LBase <= u <= LBase+LCount-1
def isCombiningV(u):
    return VBase <= VBase+VCount-1
def isCombiningT(u):
    return TBase+1 <= u <= TBase+TCount-1
def isCombinedS(u):
    return SBase <= u <= SBase+SCount-1
def decomposeS(S):
    L = (S - SBase) // NCount + LBase
    Nindex = (S - SBase) % NCount
    V = Nindex // TCount + VBase
    Tindex = Nindex % TCount
    T = Tindex + TBase if Tindex else None
    return (L,V,T)

def contourControls(contour):
    # Use second byte of the operation name (curveTo, closePath, etc),
    # as that's unique.
    return ''.join(op[0][1] for op in contour)

def outlineStructure(outline):
    return ''.join(contourControls(contour) for contour in outline)

def outlineVector(outline):
    if not outline:
        return []
    assert(outline[0][0][0] == "moveTo")
    initPos = outline[0][0][1][0]
    vec = []
    for contour in outline:
        for op in contour:
            for x_y in op[1]:
                vec.append(x_y[0] - initPos[0])
                vec.append(x_y[1] - initPos[1])
    return vec

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
    return ret

if demoS:
    S = demoS
    L,V,T = decomposeS(S)
    svgMain([fontfile, chr(L)+chr(V)+(chr(T) if T else '')+chr(S)])
    sys.exit(0)


alternates = defaultdict(list)

for weight in (100, 1000):
    shapes = {}
    split = {}
    glyphset = font.getGlyphSet(location={'wght':weight})
    for u in list(range(LBase, LBase+LCount)) + \
             list(range(VBase, VBase+VCount)) + \
             list(range(TBase+1, TBase+TCount)) + \
             list(range(SBase, SBase+SCount)):
        pen = PerContourPen(RecordingPen)
        glyphset[cmap[u]].draw(pen)
        shapes[u] = [recPen.value for recPen in pen.value]

    for S in range(SBase, SBase+SCount):
        L,V,T = decomposeS(S)
        if T is None:
            continue # Only doing LVT for now

        Llen = len(shapes[L])
        Vlen = len(shapes[V])
        Tlen = len(shapes[T])
        Slen = len(shapes[S])
        if Llen + Vlen + Tlen != Slen:
            #print("U+%04X: Contour count mismatch; skipping" % S)
            continue

        decomposed = []
        for u in (L, V, T):
            shape = shapes[u]
            component = [len(contour) for contour in shape]
            decomposed.extend(component)
        composed = [len(contour) for contour in shapes[S]]
        if decomposed != composed:
            #print("U+%04X: Contour operation count mismatch; skipping" % S)
            #continue
            pass

        # Chop shape for S into L,V,T components and save to respective lists
        # Assumption, I know...
        Sshape = shapes[S]
        Lshape = Sshape[:Llen]
        Vshape = Sshape[Llen:Llen+Vlen]
        Tshape = Sshape[Llen+Vlen:]

        alternates[L].append(Lshape)
        alternates[V].append(Vshape)
        alternates[T].append(Tshape)
        split[S] = (Lshape,Vshape,Tshape)

        #print("U+%04X: Good to go" % S)

for unicode,alts in sorted(alternates.items()):
    counter = Counter()
    sortedAlts = defaultdict(list)
    for outline in alts:
        structure = outlineStructure(outline)
        counter[structure] += 1
        sortedAlts[structure].append(outlineVector(outline))
    best = max(counter, key=lambda k: counter[k])
    print("U+%04X: Best structure matched %d out of %d instances." % (unicode, counter[best], len(alts)))

    # Build matrix for best structure
    samples = sortedAlts[best]

    # Remove duplicate samples, keeping order
    new_samples = {}
    for sample in samples:
        sample = tuple(sample)
        if sample not in new_samples:
            new_samples[sample] = 1
    samples = list(new_samples.keys())

    mat = np.matrix(samples)
    u,s,v = np.linalg.svd(mat, full_matrices=False)

    # Find number of "masters" to keep
    first = s[0] # Largest singular value
    for k in range(len(s)):
        if s[k] < first / 100:
            break

    # Truncate rank to k
    u = u[:,:k]
    s = s[:k]
    v = v[:k,:]

    reconst = np.round(u * np.diag(s) * v)
    error = reconst - mat
    print("Num masters %d max error without rounding masters %d" % (k, np.max(error)))

    # Multiply extracted features by singular values and be done with those values.
    v = np.diag(s) * v
    del s

    # v contains the list of shape-like features discovered, one in each row, and
    # u contains the combination factors of those, one row per sample.

    # Normalize range of each "axis" to 0-1; This extracts default master, and deltas.
    defaultMaster = np.zeros(np.shape(v[0]))
    for j in range(k):
        minV = np.min(u[:,j])
        maxV = np.max(u[:,j])
        diff = maxV - minV
        assert diff > 1e-5

        u[:,j] -= minV
        u[:,j] /= diff

        defaultMaster += v[j,:] * minV
        v[j,:] -= v[j,:] * minV
        v[j,:] *= diff

    defaultMaster = np.round(defaultMaster)
    deltas = np.round(v)

    # Reconstruct again, from defaultMaster+deltas
    reconst = defaultMaster + np.round(u * deltas)
    error = reconst - mat
    print("Num masters %d max error with rounding masters %d" % (k, np.max(error)))

    defaultMasterPenValues = reconstructRecordingPenValues(best, defaultMaster.tolist()[0])
    masters = [defaultMasterPenValues]
    for delta in deltas:
        values = reconstructRecordingPenValues(best, (defaultMaster+delta).tolist()[0])
        masters.append(values)

    instances = []
    for scalars in u:
        instance = np.matrix(defaultMaster)
        for scalar,delta in zip(scalars.tolist()[0],deltas):
            instance += scalar * delta
        values = reconstructRecordingPenValues(best, instance.tolist()[0])
        instances.append(values)

    SVGs = []
    for image in instances:
        rPen = RecordingPen()
        rPen.value = image
        pen = SVGPathPen(glyphset)
        rPen.replay(pen)
        commands = pen.getCommands()
        SVGs.append(commands)

    scale = .1
    with open("U+%04X.svg" % unicode, "w") as fd:

        cols = 16
        width = upem * cols
        height = upem * math.ceil(len(SVGs) / cols)

        print('<?xml version="1.0" encoding="UTF-8"?>', file=fd)
        print('<svg width="%d" height="%d" xmlns="http://www.w3.org/2000/svg">' % (width*scale, height*scale), file=fd)
        print('<rect width="100%" height="100%" fill="white"/>', file=fd)
        y = -upem * .5
        for i,commands in enumerate(SVGs):
            if i % cols == 0:
                y += upem
                x = upem * .5
            s = '<g transform="translate(%d %d) scale(%g -%g)"><path d="%s"/></g>' % (x*scale, y*scale, scale, scale, commands)
            print(s, file=fd)
            x += upem
        print('</svg>', file=fd)
