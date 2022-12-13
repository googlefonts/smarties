from contour import *
from font import *

from fontTools.ttLib import TTFont
from fontTools.misc.roundTools import otRound
from fontTools.pens.recordingPen import RecordingPen
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.svgPathPen import main as svgMain
from fontTools.pens.boundsPen import ControlBoundsPen
from fontTools.varLib.interpolatable import PerContourPen
from fontTools.ttLib.tables._g_l_y_f import Glyph
from fontTools.ttLib.tables.TupleVariation import TupleVariation
from collections import defaultdict
import struct
import math
import sys

if len(sys.argv) < 2:
    print("usage: python smarties-hangul.py NotoSerifKR-VF.otf")
    print("usage: python smarties-hangul.py NotoSansKR-VF.otf")
    sys.exit(1)

fontfile = sys.argv[1]
count = None
if len(sys.argv) > 2:
    arg = sys.argv[2]
    count = int(arg)

serif = 'serif' if fontfile.find("Serif") >= 0 else 'sans'
font = TTFont(fontfile)
upem = font['head'].unitsPerEm
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

def halve(l):
    n = len(l) // 2
    return l[:n], l[n:]

WEIGHTS = None
for axis in font['fvar'].axes:
    if axis.axisTag == 'wght':
        WEIGHTS = (axis.minValue, axis.maxValue)
        break
FAMILY_NAME = "butchered-hangul-" + serif

shapes = {}
Sbuild = {}
for weight in WEIGHTS:
    print("Font weight %d." % weight)
    Sbuild[weight] = {}
    shapes[weight] = {None: []}
    glyphset = font.getGlyphSet(location={'wght':weight})

    print("Gathering shapes.")
    for u in list(range(LBase, LBase+LCount)) + \
             list(range(VBase, VBase+VCount)) + \
             list(range(TBase+1, TBase+TCount)) + \
             list(range(SBase, SBase+SCount)):
        pen = PerContourPen(RecordingPen)
        glyphset[cmap[u]].draw(pen)
        shapes[weight][u] = [recPen.value for recPen in pen.value]


print("Gathering components.")
w0,w1 = WEIGHTS
matches = set()
alternates = defaultdict(list)
mismatch  = 0
num_matched = 0
not_matched = 0
for S in range(SBase, SBase+SCount):
    L,V,T = decomposeS(S)

    Llen = len(shapes[w0][L])
    Vlen = len(shapes[w0][V])
    Tlen = len(shapes[w0][T])
    Slen = len(shapes[w0][S])
    if Llen + Vlen + Tlen != Slen:
        mismatch += 1
        continue

    Sshape = shapes[w0][S]
    Pshape = shapes[w0][L] + shapes[w0][V] + shapes[w0][T]
    matchedOutline,_,assignment0 = matchOutline(Sshape, Pshape)

    if matchedOutline:
        pieces0 = matchedOutline[:Llen],matchedOutline[Llen:Llen+Vlen],matchedOutline[Llen+Vlen:]

        Sshape = shapes[w1][S]
        Pshape = shapes[w1][L] + shapes[w1][V] + shapes[w1][T]
        matchedOutline,_,assignment1 = matchOutline(Sshape, Pshape)
        # assert assignment0 == assignment1 # Doesn't match for Sans font. Sigh.

        pieces1 = matchedOutline[:Llen],matchedOutline[Llen:Llen+Vlen],matchedOutline[Llen+Vlen:]

        alternates[L].append((pieces0[0],pieces1[0]))
        alternates[V].append((pieces0[1],pieces1[1]))
        alternates[T].append((pieces0[2],pieces1[2]))

        num_matched += 1
        matches.add(S)

        Sbuild[S] = ((L,V,T), pieces0, pieces1)
    else:
        not_matched += 1

print("matched: %d not matched: %d contour count mismatch: %d " % (num_matched, not_matched, mismatch))
del alternates[None]


print("Learning.")
learned = {}
structs = {}
componentDefaultMaster = {}
componentMasters = {}
componentDeltas = {}
componentCoordinates = {}
for unicode,alts in sorted(alternates.items()):
    print("U+%04X: Structure matched %d." % (unicode, len(alts)))

    structure = outlineStructure(alts[0][0]) * len(WEIGHTS)
    structs[unicode] = structure
    samples = []
    for alt0,alt1 in alts:
        samples.append(outlineVector(alt0) + outlineVector(alt1))

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
        for images in data:
            for image in halve(images):
                rPen = RecordingPen()
                rPen.value = image
                pen = SVGPathPen(glyphset)
                rPen.replay(pen)
                commands = pen.getCommands()
                SVGs.append(commands)

    scale = .1
    with open("fonts/%s/svg/U+%04X.svg" % (serif, unicode), "w") as fd:

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


print("Building fonts")


style_name = "flat-original-variable"
file_name = "fonts/%s/%s-%s.ttf" % (serif,FAMILY_NAME, style_name)
print("Building %s font" % file_name)

fb = createFontBuilder(font, FAMILY_NAME, style_name, matches)

glyphSets = {}
for weight in WEIGHTS:
    glyphSets[weight] = {".notdef": Glyph()}
for S in matches:
    glyphName = cmap[S]
    pens = []
    commands = []
    for i,weight in enumerate(WEIGHTS):
        pens.append(createTTGlyphPen())
        build = Sbuild[S]
        order = build[0]
        pieces = build[i+1]
        command = []
        for piece in pieces:
            if not piece: continue
            for contour in piece:
                command.extend(contour)
        commands.append(command)
    cu2quPen = createCu2QuMultiPen(pens)
    replayCommandsThroughCu2QuMultiPen(commands, cu2quPen)
    for i,weight in enumerate(WEIGHTS):
        glyphSets[weight][glyphName] = pens[i].glyph()

glyphs, variations, axes = setupVariableFont(glyphSets, WEIGHTS)
fb.setupGlyf(glyphs)
fb.setupFvar(axes, [])
fb.setupGvar(variations)
fb.font['avar'] = font['avar']

print("Saving %s" % file_name)
fb.save(file_name)


style_name = "flat-variable"
file_name = "fonts/%s/%s-%s.ttf" % (serif,FAMILY_NAME, style_name)
print("Building %s font" % file_name)

fb = createFontBuilder(font, FAMILY_NAME, style_name, matches)

glyphSets = {}
for weight in WEIGHTS:
    glyphSets[weight] = {".notdef": Glyph()}
for S in matches:
    glyphName = cmap[S]
    pens = []
    commands = []
    for weight in WEIGHTS:
        pens.append(createTTGlyphPen())
        commands.append([])

    for unicode,piece0,piece1 in zip(*Sbuild[S]):
        if not piece0: continue
        position0 = outlinePosition(piece0)
        vector0 = outlineVector(piece0)
        position1 = outlinePosition(piece1)
        vector1 = outlineVector(piece1)

        piece01 = learned[unicode][vector0+vector1]
        piece0, piece1 = halve(piece01)

        piece0 = positionFlatOutline(piece0, position0)
        piece1 = positionFlatOutline(piece1, position1)
        commands[0].extend(piece0)
        commands[1].extend(piece1)

    cu2quPen = createCu2QuMultiPen(pens)
    replayCommandsThroughCu2QuMultiPen(commands, cu2quPen)
    for i,weight in enumerate(WEIGHTS):
        glyphSets[weight][glyphName] = pens[i].glyph()

glyphs, variations, axes = setupVariableFont(glyphSets, WEIGHTS)
fb.setupGlyf(glyphs)
fb.setupFvar(axes, [])
fb.setupGvar(variations)
fb.font['avar'] = font['avar']

print("Saving %s" % file_name)
fb.save(file_name)


style_name = "smarties-variable"
file_name = "fonts/%s/%s-%s.ttf" % (serif,FAMILY_NAME, style_name)
print("Building %s font" % file_name)

components = []
componentNames = {}
for unicode in learned.keys():
    # Give name to each learned item:
    name = "uni%04X" % unicode
    componentNames[unicode] = name
    components.append(name)

fb = createFontBuilder(font, FAMILY_NAME, style_name, matches, components)

variations = {}

# Write out components & gather variations
glyphs = {".notdef": Glyph()}
for unicode in learned.keys():
    glyphName = componentNames[unicode]
    deltas = componentDeltas[unicode]
    variations[glyphName] = []

    masterCommands = componentMasters[unicode]
    # split masterCommands into their weight components
    master0s = []
    master1s = []
    for master in masterCommands:
        master0,master1 = halve(master)
        master0s.append(master0)
        master1s.append(master1)

    numMasters = len(masterCommands)
    pens = [createTTGlyphPen() for i in range(numMasters * len(WEIGHTS))]
    # Replay all masters together through cu2qu multi-pen!
    cu2quPen = createCu2QuMultiPen(pens)
    replayCommandsThroughCu2QuMultiPen(master0s + master1s, cu2quPen)

    masterGlyph = pens[0].glyph()
    weightGlyph = pens[numMasters].glyph()
    glyphs[glyphName] = masterGlyph
    allDeltas = []

    for i,pen in enumerate(pens):
        if i == 0: # Skip default master
            allDeltas.append(None)
            continue
        # Weight master is already loaded; we can't load it again
        glyph = pen.glyph() if i != numMasters else weightGlyph
        # Subtract base; either the default master or the weight master
        coords = glyph.coordinates - (masterGlyph.coordinates if i <= numMasters else weightGlyph.coordinates)
        allDeltas.append(coords.copy())
        if i > numMasters:
            coords -= allDeltas[i - numMasters]
        coords.extend([(0,0), (0,0), (0,0), (0,0)]) # TODO Phantom points
        axes = {}
        if i % numMasters != 0:
            tag = "%04d" % ((i % numMasters) - 1)
            axes[tag] = (0, 1, 1)
        if i >= numMasters:
            axes['wght'] = (0, 1, 1)

        tv = TupleVariation(axes, coords)
        variations[glyphName].append(tv)


# Write out composites.
reverseGlyphMap = fb.font.getReverseGlyphMap()
for S in matches:
    glyphName = cmap[S]
    order0,pieces0,pieces1 = Sbuild[S]
    glyph = Glyph()

    boundsPen = ControlBoundsPen(None)
    rPen = RecordingPen()
    for order,piece in zip(order0,pieces0):
        for contour in piece:
            rPen.value.extend(contour)
    rPen.replay(boundsPen)
    b = [otRound(v) for v in boundsPen.bounds]
    data = bytearray(struct.pack(">hhhhh", -2, b[0], b[1], b[2], b[3]))
    variation = []
    for componentUnicode,piece0,piece1 in zip(order0,pieces0,pieces1):
        if not piece0: continue
        componentName = componentNames[componentUnicode]

        position0 = outlinePosition(piece0)
        position0 = [otRound(v) for v in position0]
        position1 = outlinePosition(piece1)
        position1 = [otRound(v) for v in position1]

        vector0 = outlineVector(piece0)
        vector1 = outlineVector(piece1)
        piece = learned[componentUnicode][vector0+vector1]
        vector = outlineVector(piece, flat=True)
        coordinates = componentCoordinates[componentUnicode][vector]

        # Build glyph data

        flag = struct.pack(">H", (1<<3)|(1<<4))
        numAxes = struct.pack(">B", len(coordinates))
        gid = struct.pack(">H", reverseGlyphMap[componentName])
        axisIndices = b''.join(struct.pack(">B", i) for i in range(len(coordinates)))
        axisValues = b''.join(struct.pack(">H", otRound(v * 16384)) for v in coordinates)
        translate = struct.pack(">hh", *position0)

        rec = flag + numAxes + gid + axisIndices + axisValues + translate
        data.extend(rec)

        # Build variation

        x,y = position1[0] - position0[0], position1[1] - position0[1]
        variation.append((x, y)) # Translate

    glyph.data = bytes(data)
    glyphs[glyphName] = glyph

    variation.extend([(0,0), (0,0), (0,0), (0,0)]) # Phantom points TODO
    axes = {'wght': (0, 1, 1)}
    tv = TupleVariation(axes, variation)
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
    fb.font['avar'].segments[axis.axisTag] = {}

fb.setupGvar(variations)

fb.font.recalcBBoxes = False

print("Saving %s" % file_name)
fb.save(file_name)
