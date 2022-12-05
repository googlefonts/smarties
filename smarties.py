
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import RecordingPen
from fontTools.pens.svgPathPen import main as svgMain
from fontTools.varLib.interpolatable import PerContourPen
from collections import Counter, defaultdict
import numpy as np
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
cmap = font['cmap'].getBestCmap()
glyphset = font.getGlyphSet()

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
    return tuple(contourControls(contour) for contour in outline)

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

if demoS:
    S = demoS
    L,V,T = decomposeS(S)
    svgMain([fontfile, chr(L)+chr(V)+(chr(T) if T else '')+chr(S)])
    sys.exit(0)


shapes = {}
alternates = defaultdict(list)
split = {}

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

for u,alts in sorted(alternates.items()):
    counter = Counter()
    sortedAlts = defaultdict(list)
    for outline in alts:
        structure = outlineStructure(outline)
        counter[structure] += 1
        sortedAlts[structure].append(outlineVector(outline))
    best = max(counter, key=lambda k: counter[k])
    print("U+%04X: Best structure matched %d out of %d instances." % (u, counter[best], len(alts)))

    # Build matrix for best structure
    samples = sortedAlts[best]
    mat = np.transpose(np.matrix(samples))
    ret = np.linalg.svd(mat, full_matrices=False)
    print(ret)
    break
