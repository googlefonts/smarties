
from fontTools.ttLib import TTFont
from fontTools.misc.roundTools import otRound
from fontTools.pens.recordingPen import RecordingPen
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.svgPathPen import main as svgMain
from fontTools.pens.statisticsPen import StatisticsPen
from fontTools.varLib.interpolatable import PerContourPen
from fontTools.fontBuilder import FontBuilder
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.pens.cu2quPen import Cu2QuPen, Cu2QuMultiPen
from fontTools.ttLib.tables._g_l_y_f import Glyph, GlyphCoordinates, GlyphComponent
from fontTools.ttLib.tables.TupleVariation import TupleVariation
from collections import defaultdict
from itertools import permutations
import numpy as np
from scipy.optimize import linear_sum_assignment
import struct
import math
import sys

if len(sys.argv) < 2:
    print("usage: python smarties.py SourceHanSerifKR-VF.otf")
    sys.exit(1)

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

def decomposeS(S):
    L = (S - SBase) // NCount + LBase
    Nindex = (S - SBase) % NCount
    V = Nindex // TCount + VBase
    Tindex = Nindex % TCount
    T = Tindex + TBase if Tindex else None
    return (L,V,T)

def contourStructure(contour):
    # Use second byte of the operation name (curveTo, closePath, etc),
    # as that's unique.
    return ''.join(op[0][1] for op in contour)

def outlineStructure(outline):
    return ''.join(contourStructure(contour) for contour in outline)

def outlinePosition(outline):
    assert(outline[0][0][0] == "moveTo")
    return outline[0][0][1][0]

def flatOutlinePosition(outline, initPos):
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

def matchOutline(shape, ref):
    assert len(shape) == len(ref)

    # Shortcut: If structures match assume it's correct.
    # Although if order is wrong matching can recover the
    # correct order...
    #if outlineStructure(shape) == outlineStructure(ref):
    #    return shape, 0

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
        return None, 1e10

    # We have a matching. Reorder contours and return
    reordered = []
    for c in cols:
        reordered.append(shape[c])

    return reordered, cost



if demoS:
    S = demoS
    L,V,T = decomposeS(S)
    svgMain([fontfile, chr(L)+chr(V)+(chr(T) if T else '')+chr(S)])
    sys.exit(0)


alternates = defaultdict(list)
matches = set()
Sbuild = {}
WEIGHTS = (250, 900)

for weight in WEIGHTS:
    print("Font weight %d." % weight)
    Sbuild[weight] = {}
    mismatch  = 0
    num_matched = 0
    not_matched = 0
    shapes = {}
    glyphset = font.getGlyphSet(location={'wght':weight})

    print("Gathering shapes.")
    for u in list(range(LBase, LBase+LCount)) + \
             list(range(VBase, VBase+VCount)) + \
             list(range(TBase+1, TBase+TCount)) + \
             list(range(SBase, SBase+SCount)):
        pen = PerContourPen(RecordingPen)
        glyphset[cmap[u]].draw(pen)
        shapes[u] = [recPen.value for recPen in pen.value]

    print("Gathering components.")
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
            mismatch += 1
            continue

        Sshape = shapes[S]
        bestCost = 1e10
        bestOrder = None
        bestOutlines = None
        for order in permutations((L,V,T)):
            # Chop shape for S into L,V,T components and save to respective lists
            # Assumption, I know...
            len0 = len(shapes[order[0]])
            len1 = len(shapes[order[1]])
            len2 = len(shapes[order[2]])
            shape0 = Sshape[:len0]
            shape1 = Sshape[len0:len0+len1]
            shape2 = Sshape[len0+len1:]

            matchedOutline0,cost0 = matchOutline(shape0, shapes[order[0]])
            matchedOutline1,cost1 = matchOutline(shape1, shapes[order[1]])
            matchedOutline2,cost2 = matchOutline(shape2, shapes[order[2]])
            if cost0 + cost1 + cost2 < bestCost:
               bestOrder = order
               bestCost = cost0 + cost1 + cost2
               bestOutlines = matchedOutline0,matchedOutline1,matchedOutline2

        if bestOrder:
            alternates[bestOrder[0]].append(bestOutlines[0])
            alternates[bestOrder[1]].append(bestOutlines[1])
            alternates[bestOrder[2]].append(bestOutlines[2])
            num_matched += 1
            matches.add(S)
            Sbuild[weight][S] = (bestOrder, (bestOutlines[0], bestOutlines[1], bestOutlines[2]))
        else:
            not_matched += 1

    print("matched: %d not matched: %d mismatch: %d " % (num_matched, not_matched, mismatch))

print("Learning.")
learned = {}
structs = {}
componentDefaultMaster = {}
componentMasters = {}
componentDeltas = {}
componentCoordinates = {}
for unicode,alts in sorted(alternates.items()):
    print("U+%04X: Structure matched %d." % (unicode, len(alts)))

    structure = outlineStructure(alts[0])
    structs[unicode] = structure
    samples = []
    for alt in alts:
        samples.append(outlineVector(alt))

    # Remove duplicate samples, keeping order
    new_samples = {}
    for sample in samples:
        if sample not in new_samples:
            new_samples[sample] = 1
    samples = list(new_samples.keys())

    mat = np.matrix(samples)
    u,s,v = np.linalg.svd(mat, full_matrices=False)

    # Find number of "masters" to keep
    first = s[0] # Largest singular value
    k = len(s)
    while k and s[k - 1] < first / 100:
        k -= 1

    # Truncate rank to k
    u = u[:,:k]
    s = s[:k]
    v = v[:k,:]

    reconst = np.round(u * np.diag(s) * v)
    error = reconst - mat
    maxError = np.max(error)
    meanSqError = np.mean(np.square(error))
    #print("Num masters %d max error %d mean-squared error %g" % (k, maxError, meanSqError))

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

        u[:,j] -= minV
        if diff:
            u[:,j] /= diff

        defaultMaster += v[j,:] * minV
        v[j,:] *= diff

    defaultMaster = np.round(defaultMaster)
    deltas = np.round(v)
    # Round scalars to 2.14
    u = np.round(u * 16384) / 16384

    # Reconstruct again, from defaultMaster+deltas
    reconst = defaultMaster + u * deltas
    error = reconst - mat
    maxError = np.max(error)
    meanSqError = np.mean(np.square(error))
    print("Num masters %d max error %d mean-squared error %g" % (k+1, maxError, meanSqError))

    defaultMasterPenValues = reconstructRecordingPenValues(structure, defaultMaster.tolist()[0])
    componentDefaultMaster[unicode] = defaultMasterPenValues
    masters = [defaultMasterPenValues]
    componentDeltas[unicode] = []
    for delta in deltas:
        componentDeltas[unicode].append(reconstructRecordingPenValues(structure, delta.tolist()[0]))
        values = reconstructRecordingPenValues(structure, (defaultMaster+delta).tolist()[0])
        masters.append(values)
    componentMasters[unicode] = masters


    instances = []
    componentCoordinates[unicode] = {}
    for scalars in u:
        instance = np.matrix(defaultMaster)
        scals = scalars.tolist()[0]
        for scalar,delta in zip(scals,deltas):
            instance += scalar * delta
        instance = np.round(instance)
        instance = tuple(instance.tolist()[0])
        componentCoordinates[unicode][instance] = scals
        values = reconstructRecordingPenValues(structure, instance)
        instances.append(values)

    learned[unicode] = {}
    for s,i in zip(samples,instances):
        learned[unicode][s] = i


    unique_instances = set(tuple(i) for i in instances)
    print("Num instances %d num unique instances %d" % (len(instances), len(unique_instances)))
    del unique_instances

    originals = []
    for sample in samples:
        values = reconstructRecordingPenValues(structure, sample)
        originals.append(values)

    masterSVGs = []
    instanceSVGs = []
    origSVGs = []
    for data,SVGs in ((masters,masterSVGs), (instances,instanceSVGs), (originals,origSVGs)):
        for image in data:
            rPen = RecordingPen()
            rPen.value = image
            pen = SVGPathPen(glyphset)
            rPen.replay(pen)
            commands = pen.getCommands()
            SVGs.append(commands)

    scale = .1
    with open("U+%04X.svg" % unicode, "w") as fd:

        cols = 16
        width = upem * (cols + 1)
        height = upem * (math.ceil(len(masterSVGs) / cols) + math.ceil(len(instanceSVGs) / cols) + 1)

        print('<?xml version="1.0" encoding="UTF-8"?>', file=fd)
        print('<svg width="%d" height="%d" xmlns="http://www.w3.org/2000/svg">' % (width*scale, height*scale), file=fd)
        print('<rect width="100%" height="100%" fill="white"/>', file=fd)
        y = 0
        for i,commands in enumerate(masterSVGs):
            if i % cols == 0:
                y += upem
                x = upem
            s = '<g fill="green" transform="translate(%d %d) scale(%g -%g)"><path d="%s"/></g>' % (x*scale, y*scale, scale, scale, commands)
            print(s, file=fd)
            x += upem
        for i,(origCommands,instCommands) in enumerate(zip(origSVGs,instanceSVGs)):
            if i % cols == 0:
                y += upem
                x = upem
            s = '<g fill="red" transform="translate(%d %d) scale(%g -%g)"><path d="%s"/></g>' % (x*scale, y*scale, scale, scale, origCommands)
            s += '\n'
            s += '<g fill="black" opacity=".8" transform="translate(%d %d) scale(%g -%g)"><path d="%s"/></g>' % (x*scale, y*scale, scale, scale, instCommands)
            print(s, file=fd)
            x += upem
        print('</svg>', file=fd)

def createFontBuilder(font, style, chars, extraGlyphs=[]):
    upem = font['head'].unitsPerEm
    cmap = font['cmap'].getBestCmap()
    subset_cmap = {u:g for u,g in cmap.items() if u in chars}
    subset_glyphs = set(subset_cmap.values())
    subset_glyphOrder = [g for g in font.getGlyphOrder() if g in subset_glyphs] + extraGlyphs
    subset_glyphOrder.insert(0, ".notdef")
    del subset_glyphs
    metrics = font['hmtx'].metrics
    subset_metrics = {g:metrics[g] if g in metrics else (0,0) for g in subset_glyphOrder}
    nameStrings = dict(
        familyName=dict(en="butchered-hangul-serif"),
        styleName=dict(en=style),
    )

    fb = FontBuilder(upem, isTTF=True)
    fb.setupNameTable(nameStrings)
    fb.setupGlyphOrder(subset_glyphOrder)
    fb.setupCharacterMap(subset_cmap)
    fb.setupHorizontalMetrics(subset_metrics)
    hhea = font['hhea']
    fb.setupHorizontalHeader(ascent=hhea.ascent, descent=hhea.descent)
    os2 = font['OS/2']
    fb.setupOS2(sTypoAscender=os2.sTypoAscender, usWinAscent=os2.usWinAscent, usWinDescent=os2.usWinDescent)
    fb.setupPost()

    return fb

def createCu2QuPen(pen):
    return Cu2QuPen(pen, .5, reverse_direction=True)
def createCu2QuMultiPen(pens):
    return Cu2QuMultiPen(pens, .5)#, reverse_direction=True) # https://github.com/fonttools/fonttools/issues/2914

def replayCommandsThroughCu2QuMultiPen(commands, cu2quPen):
    for ops in zip(*commands):
        opNames = [op[0] for op in ops]
        opArgs = [op[1] for op in ops]
        opName = opNames[0]
        assert all(name == opName for name in opNames)
        if len(opArgs[0]):
            getattr(cu2quPen, opName)(opArgs)
        else:
            getattr(cu2quPen, opName)()

def getCoords(glyph):
    if glyph.numberOfContours > 0:
        return glyph.coordinates
    if glyph.numberOfContours == 0:
        return GlyphCoordinates()
    if glyph.numberOfContours == -1:
        components = glyph.components
        coords = []
        for comp in components:
            coords.append((comp.x, comp.y))
        return GlyphCoordinates(coords)

    assert False

def setupVariableFont(glyphSets):
    variations = {}
    tag = 'wght'
    glyphs = glyphSets[WEIGHTS[0]]
    varGlyphs = glyphSets[WEIGHTS[1]]
    for glyphName in glyphs:
        glyph = glyphs[glyphName]
        varGlyph = varGlyphs[glyphName]

        coords = getCoords(varGlyph) - getCoords(glyph)

        # Add phantom points
        coords.extend([(0,0), (0,0), (0,0), (0,0)])
        axes = {tag: (0, 1, 1)}
        tv = TupleVariation(axes, coords)
        variations[glyphName] = [tv]

    axes = [(tag, WEIGHTS[0], WEIGHTS[0], WEIGHTS[1], tag)]
    return glyphs, variations, axes



print("Building fonts")


print("Building butchered-hangul-serif-flat-original-variable font")

fb = createFontBuilder(font, "flat-original-variable", matches)
glyphSets = {}
for weight in WEIGHTS:
    glyphSets[weight] = {".notdef": Glyph()}
for S in matches:
    glyphName = cmap[S]
    pens = []
    commands = []
    for weight in WEIGHTS:
        pens.append(TTGlyphPen(None))
        order,pieces = Sbuild[weight][S]
        command = []
        for piece in pieces:
            for contour in piece:
                command.extend(contour)
        commands.append(command)
    cu2quPen = createCu2QuMultiPen(pens)
    replayCommandsThroughCu2QuMultiPen(commands, cu2quPen)
    for i,weight in enumerate(WEIGHTS):
        glyphSets[weight][glyphName] = pens[i].glyph()

glyphs, variations, axes = setupVariableFont(glyphSets)
fb.setupGlyf(glyphs)
fb.setupFvar(axes, [])
fb.setupGvar(variations)
fb.font['avar'] = font['avar']

print("Saving butchered-hangul-serif-flat-original-variable.ttf")
fb.save("butchered-hangul-serif-flat-original-variable.ttf")


print("Building butchered-hangul-serif-flat-variable font")
fb = createFontBuilder(font, "flat-variable", matches)
glyphSets = {}
for weight in WEIGHTS:
    glyphSets[weight] = {".notdef": Glyph()}
for S in matches:
    glyphName = cmap[S]
    pens = []
    commands = []
    for weight in WEIGHTS:
        pens.append(TTGlyphPen(None))
        order,pieces = Sbuild[weight][S]
        command = []
        for unicode,piece in zip(order,pieces):
            position = outlinePosition(piece)
            vector = outlineVector(piece)
            piece = learned[unicode][vector]
            piece = flatOutlinePosition(piece, position)
            command.extend(piece)
        commands.append(command)
    cu2quPen = createCu2QuMultiPen(pens)
    replayCommandsThroughCu2QuMultiPen(commands, cu2quPen)
    for i,weight in enumerate(WEIGHTS):
        glyphSets[weight][glyphName] = pens[i].glyph()

glyphs, variations, axes = setupVariableFont(glyphSets)
fb.setupGlyf(glyphs)
fb.setupFvar(axes, [])
fb.setupGvar(variations)
fb.font['avar'] = font['avar']

print("Saving butchered-hangul-serif-flat-variable.ttf")
fb.save("butchered-hangul-serif-flat-variable.ttf")


print("Building butchered-hangul-serif-smarties-variable font")
components = []
componentNames = {}
for unicode in learned.keys():
    # Give name to each learned item:
    name = "uni%04X" % unicode
    componentNames[unicode] = name
    components.append(name)

fb = createFontBuilder(font, "smarties-variable", matches, components)
variations = {}

# Write out components & gather variations
glyphs = {".notdef": Glyph()}
for unicode in learned.keys():
    glyphName = componentNames[unicode]
    deltas = componentDeltas[unicode]
    variations[glyphName] = []

    masterCommands = componentMasters[unicode]
    numMasters = len(masterCommands)
    pens = [TTGlyphPen(None) for i in range(numMasters)]
    # Replay all masters together through cu2qu multi-pen!
    cu2quPen = createCu2QuMultiPen(pens)
    replayCommandsThroughCu2QuMultiPen(masterCommands, cu2quPen)

    masterGlyph = pens[0].glyph()
    glyphs[glyphName] = masterGlyph
    for i,pen in enumerate(pens[1:]):
        glyph = pen.glyph()
        coords = glyph.coordinates - masterGlyph.coordinates
        # Add phantom points
        coords.extend([(0,0), (0,0), (0,0), (0,0)])
        tag = "%04d" % i
        axes = {tag: (0, 1, 1)}

        tv = TupleVariation(axes, coords)
        variations[glyphName].append(tv)

# Write out composites.
reverseGlyphMap = fb.font.getReverseGlyphMap()
for S in matches:
    glyphName = cmap[S]
    order0,pieces0 = Sbuild[WEIGHTS[0]][S]
    order1,pieces1 = Sbuild[WEIGHTS[1]][S]
    assert order0 == order1
    glyph = Glyph()

    data = bytearray(struct.pack(">h", -2) + b"\00\00\00\00\00\00\00\00")
    variation = []
    for componentUnicode,piece0,piece1 in zip(order0,pieces0,pieces1):
        componentName = componentNames[componentUnicode]

        position0 = outlinePosition(piece0)
        position0 = [otRound(v) for v in position0]
        position1 = outlinePosition(piece1)
        position1 = [otRound(v) for v in position1]

        vector = outlineVector(piece0)
        piece0 = learned[componentUnicode][vector]
        vector = outlineVector(piece0, flat=True)
        coordinates0 = componentCoordinates[componentUnicode][vector]
        vector = outlineVector(piece1)
        piece1 = learned[componentUnicode][vector]
        vector = outlineVector(piece1, flat=True)
        coordinates1 = componentCoordinates[componentUnicode][vector]
        assert len(coordinates0) == len(coordinates1)

        # Build glyph data

        flag = struct.pack(">H", (1<<3)|(1<<4))
        gid = struct.pack(">H", reverseGlyphMap[componentName])
        numAxes = struct.pack(">B", len(coordinates0))
        axisIndices = b''.join(struct.pack(">B", i) for i in range(len(coordinates0)))
        axisValues = b''.join(struct.pack(">H", otRound(v * 16384)) for v in coordinates0)
        translate = struct.pack(">hh", *position0)

        rec = flag + gid + numAxes + axisIndices + axisValues + translate
        data.extend(rec)

        # Build variation

        for coord0,coord1 in zip(coordinates0, coordinates1):
            delta = otRound(coord1 * 16384) - otRound(coord0 * 16384)
            variation.append((delta, 0))
        x,y = position1[0] - position0[0], position1[1] - position0[1]
        variation.append((x, y)) # Translate
        variation.append((0, 0)) # Rotation
        variation.append((0, 0)) # Scale
        variation.append((0, 0)) # Skew
        variation.append((0, 0)) # TCenter

    glyph.data = bytes(data)
    glyphs[glyphName] = glyph

    variation.extend([(0,0), (0,0), (0,0), (0,0)]) # Phantom points
    axes = {'wght': (0, 1, 1)}
    tv = TupleVariation(axes, coords)
    variations[glyphName] = [tv]

fb.setupGlyf(glyphs)

# Setup fvar/gvar
numAxes = 0
for coordArray in componentCoordinates.values():
    numAxes = max(numAxes, len(next(iter(coordArray.values()))))
axes = []
for i in range(numAxes):
    tag = "%04d" % i
    axes.append((tag, 0, 0, 1, tag))
axes.append(('wght', WEIGHTS[0], WEIGHTS[0], WEIGHTS[1], "Weight"))
fb.setupFvar(axes, [])
fb.font['avar'] = font['avar']
# Hide axes and add avar segments
for axis in fb.font['fvar'].axes:
    if axis.axisTag == 'wght': continue
    axis.flags = 1 # HIDDEN_AXIS
    fb.font['avar'].segments[axis.axisTag] = {-1.0: -1.0, 0.0: 0.0, 1.0: 1.0}

fb.setupGvar(variations)

print("Saving butchered-hangul-serif-smarties-variable.ttf")
fb.font.recalcBBoxes = False
fb.save("butchered-hangul-serif-smarties-variable.ttf", )
