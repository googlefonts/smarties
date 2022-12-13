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
from itertools import product
import re
import struct
import math
import sys

if len(sys.argv) < 2:
    print("usage: python smarties-hangul.py NotoSerifCJKtc-VF.otf")
    print("usage: python smarties-hangul.py NotoSansCJKtc-VF.otf")
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

# Read database
kangxiUnicodes = {}
Hbuild = {}
with open("Unihan_RadicalStrokeCounts.txt") as f:
    for line in f:
        if count is not None:
            count -= 1
            if not count:
                break
        if line.startswith("#"): continue
        fields = line.split()
        if not fields: continue
        if fields[1] != "kRSAdobe_Japan1_6": continue
        unicode = int(fields[0][2:], 16)
        if unicode not in cmap: continue
        radicals = [f for f in fields[2:] if f[0] == 'C']
        Hbuild[unicode] = []
        # Parse radicals
        cid = None
        for radical in radicals:
            matches = re.match("C\+([0-9]*)\+([0-9]*)\.([0-9]*)\.([0-9]*)", radical)

            if cid is None:
                cid = matches[1]
            else:
                assert cid == matches[1], (cid, matches[1])

            kangxi = int(matches[2])
            kangxiStrokes = int(matches[3])
            otherStrokes = int(matches[4])
            Hbuild[unicode].append((kangxi, kangxiStrokes, otherStrokes))

            if len(radicals) == 1 and otherStrokes == 0:
                if kangxi not in kangxiUnicodes:
                    kangxiUnicodes[kangxi] = []
                kangxiUnicodes[kangxi].append(unicode)

WEIGHTS = None
for axis in font['fvar'].axes:
    if axis.axisTag == 'wght':
        WEIGHTS = (axis.minValue, axis.maxValue)
        break
FAMILY_NAME = "butchered-han-" + serif

shapes = {}
for weight in WEIGHTS:
    print("Font weight %d." % weight)
    shapes[weight] = {None: []}
    glyphset = font.getGlyphSet(location={'wght':weight})

    print("Gathering shapes.")
    for u in Hbuild:
        pen = PerContourPen(RecordingPen)
        glyphset[cmap[u]].draw(pen)
        shapes[weight][u] = [recPen.value for recPen in pen.value]

print("Gathering components.")
w0,w1 = WEIGHTS
matches = set()
num_matched = 0
not_matched = 0
alternates = defaultdict(list)
otherAlternates = defaultdict(list)
Hbuild2 = {}
for H in Hbuild:
    Hshape0 = shapes[w0][H]
    Hshape1 = shapes[w0][H]
    kangxis = [b[0] for b in Hbuild[H]]
    kangxiUnicodeLists = []
    for kangxi in kangxis:
        kangxiUnicodeLists.append(kangxiUnicodes.get(kangxi, []))

    bestMatchedOutline0 = None
    bestCost = 1e20
    bestAssignment = None
    bestUnicodeList = None
    for unicodeList in product(*kangxiUnicodeLists):
        # Collect radical shapes
        Rshape = []
        for radical in unicodeList:
            Rshape.extend(shapes[w0][radical])

        if len(Hshape0) < len(Rshape): continue

        matchedOutline,cost,assignment = matchOutline(Hshape0, Rshape, partial=True)

        if matchedOutline and cost < bestCost:
            bestMatchedOutline0 = matchedOutline
            bestCost = cost
            bestAssignment = assignment
            bestUnicodeList = unicodeList

    if bestMatchedOutline0:
        num_matched += 1
        matches.add(H)

        bestMatchedOutline1 = reorderAssignment(shapes[w1][H], bestAssignment)
        assert len(bestMatchedOutline0) == len(bestMatchedOutline1)

        offset = 0;
        rShapes = []
        for radical in bestUnicodeList:
            baseRshape = shapes[w0][radical]
            Rlen = len(baseRshape)
            rShape0 = bestMatchedOutline0[offset:offset+Rlen]
            rShape1 = bestMatchedOutline1[offset:offset+Rlen]
            offset += Rlen

            rShapes.append((radical, rShape0, rShape1))
            alternates[radical].append((rShape0, rShape1))

        otherStrokes0 = bestMatchedOutline0[offset:]
        otherStrokes1 = bestMatchedOutline1[offset:]
        structure = ''
        if otherStrokes0:
            structure = outlineStructure(otherStrokes0)
            otherAlternates[structure].append((otherStrokes0, otherStrokes1))

            rShapes.append((structure, otherStrokes0, otherStrokes1))

        Hbuild2[H] = rShapes

    else:
        not_matched += 1

print("matched: %d not matched: %d." % (num_matched, not_matched))

print("%d kangxi matched %d instances." % (len(alternates), sum(len(a) for a in alternates.values())))
print("%d other-strokes." % len(otherAlternates))

alternates.update(otherAlternates)
del otherAlternates

print("Learning.")
learned = {}
structs = {}
componentDefaultMaster = {}
componentMasters = {}
componentDeltas = {}
componentCoordinates = {}
for key,alts in alternates.items():
    if type(key) == int:
        print("U+%04X: Structure matched %d." % (key, len(alts)))

    structure = outlineStructure(alts[0][0]) * len(WEIGHTS)
    structs[key] = structure
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
    componentDefaultMaster[key] = defaultMasterPenValues
    masters = [defaultMasterPenValues]
    componentDeltas[key] = []
    for delta in deltas:
        componentDeltas[key].append(reconstructRecordingPenValues(structure, delta.tolist()[0]))
        values = reconstructRecordingPenValues(structure, (defaultMaster+delta).tolist()[0])
        masters.append(values)
    componentMasters[key] = masters

    instances = []
    componentCoordinates[key] = {}
    for scalars in u:
        instance = np.matrix(defaultMaster)
        scals = scalars.tolist()[0]
        for scalar,delta in zip(scals,deltas):
            instance += scalar * delta
        instance = np.round(instance)
        instance = tuple(instance.tolist()[0])
        componentCoordinates[key][instance] = scals
        values = reconstructRecordingPenValues(structure, instance)
        instances.append(values)

    learned[key] = {}
    for s,i in zip(samples,instances):
        learned[key][s] = i


    unique_instances = set(tuple(i) for i in instances)
    print("Num instances %d num unique instances %d" % (len(instances), len(unique_instances)))
    del unique_instances

    if type(key) != int:
        continue

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
    with open("fonts/%s/svg/han/U+%04X.svg" % (serif, key), "w") as fd:

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
print("Building %s" % file_name)

fb = createFontBuilder(font, FAMILY_NAME, style_name, matches)

glyphSets = {}
for weight in WEIGHTS:
    glyphSets[weight] = {".notdef": Glyph()}
for H in matches:
    glyphName = cmap[H]
    pens = []
    commands = []
    for i,weight in enumerate(WEIGHTS):
        pens.append(createTTGlyphPen())
        rShapes = Hbuild2[H]

        command = []
        for rShape in rShapes:
            shape = rShape[i+1]
            for contour in shape:
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
print("Building %s" % file_name)

fb = createFontBuilder(font, FAMILY_NAME, style_name, matches)

glyphSets = {}
for weight in WEIGHTS:
    glyphSets[weight] = {".notdef": Glyph()}
for H in matches:
    glyphName = cmap[H]
    pens = []
    commands = []
    for i,weight in enumerate(WEIGHTS):
        pens.append(createTTGlyphPen())
        commands.append([])

    """
        command = []
        for rShape in rShapes:
            shape = rShape[i+1]

            position1 = outlinePosition(piece1)
            vector1 = outlineVector(piece1)

            piece01 = learned[unicode][vector0+vector1]
            piece0, piece1 = halve(piece01)

            piece0 = positionFlatOutline(piece0, position0)
            piece1 = positionFlatOutline(piece1, position1)
            commands[0].extend(piece0)

            for contour in shape:
                command.extend(contour)
    """

    for key,piece0,piece1 in Hbuild2[H]:
        position0 = outlinePosition(piece0)
        vector0 = outlineVector(piece0)
        position1 = outlinePosition(piece1)
        vector1 = outlineVector(piece1)

        piece01 = learned[key][vector0+vector1]
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
print("Building %s" % file_name)

components = []
componentNames = {}
i = 0
for key in learned.keys():
    # Give name to each learned item:
    if type(key) == int:
        name = "uni%04X" % key
    else:
        name = "comp%d" % i
        i += 1
    componentNames[key] = name
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
for H in matches:
    glyphName = cmap[H]
    keys = []
    pieces0 = []
    pieces1 = []
    for key,piece0,piece1 in Hbuild2[H]:
        pieces0.append(piece0)
        pieces1.append(piece1)

    glyph = Glyph()

    boundsPen = ControlBoundsPen(None)
    rPen = RecordingPen()
    for _,piece,_ in Hbuild2[H]:
        for contour in piece:
            rPen.value.extend(contour)
    rPen.replay(boundsPen)
    b = [otRound(v) for v in boundsPen.bounds]
    data = bytearray(struct.pack(">hhhhh", -2, b[0], b[1], b[2], b[3]))
    variation = []
    for componentKey,piece0,piece1 in Hbuild2[H]:
        if not piece0: continue
        componentName = componentNames[componentKey]

        position0 = outlinePosition(piece0)
        position0 = [otRound(v) for v in position0]
        position1 = outlinePosition(piece1)
        position1 = [otRound(v) for v in position1]

        vector0 = outlineVector(piece0)
        vector1 = outlineVector(piece1)
        piece = learned[componentKey][vector0+vector1]
        vector = outlineVector(piece, flat=True)
        coordinates = componentCoordinates[componentKey][vector]

        # Work around our 1-master 2-master issue.
        if coordinates == [0]:
            coordinates = []

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
