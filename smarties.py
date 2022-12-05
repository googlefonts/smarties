
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import RecordingPen
from fontTools.pens.svgPathPen import main as svgMain
from fontTools.varLib.interpolatable import PerContourPen
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


if demoS:
    S = demoS
    L,V,T = decomposeS(S)
    svgMain([fontfile, chr(L)+chr(V)+(chr(T) if T else '')+chr(S)])
    sys.exit(0)


shapes = {}

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

    if len(shapes[L])+len(shapes[V])+len(shapes[T]) != len(shapes[S]):
        print("U+%04X: Contour count mismatch; skipping" % S)
        continue

    decomposed = []
    for u in (L, V, T):
        shape = shapes[u]
        component = [len(contour) for contour in shape]
        print(component)
        decomposed.extend(component)
    composed = [len(contour) for contour in shapes[S]]
    print(decomposed)
    print(composed)
    if decomposed != composed:
        print("U+%04X: Contour operation count mismatch; skipping" % S)
        continue

    print("U+%04X: Good to go" % S)
