from contour import *
from font import *

from fontTools.ttLib import TTFont
from fontTools.misc.roundTools import otRound
from fontTools.misc.transform import Transform, Identity
from fontTools.pens.recordingPen import RecordingPen
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.svgPathPen import main as svgMain
from fontTools.pens.boundsPen import ControlBoundsPen
from fontTools.varLib.interpolatable import PerContourPen
from fontTools.ttLib.tables._g_l_y_f import Glyph
from fontTools.ttLib.tables.TupleVariation import TupleVariation
from collections import defaultdict
from enum import IntEnum
import re
import struct
import math
import sys

if len(sys.argv) < 2:
    print("usage: python smarties-han.py NotoSerifCJKtc-VF.otf")
    print("usage: python smarties-han.py NotoSansCJKtc-VF.otf")
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

def isCircledNumber(u):
    if type(u) == str:
        u = ord(u)
    return 0x2460 <= u <= 0x2473


class IdeoDescription(IntEnum):
    LEFT_TO_RIGHT               = 0x2FF0 # ⿰
    ABOVE_TO_BELOW              = 0x2FF1 # ⿱
    LEFT_TO_MIDDLE_AND_RIGHT    = 0x2FF2 # ⿲
    ABOVE_TO_MIDDLE_AND_BELOW   = 0x2FF3 # ⿳
    FULL_SURROUND               = 0x2FF4 # ⿴
    SURROUND_FROM_ABOVE         = 0x2FF5 # ⿵
    SURROUND_FROM_BELOW         = 0x2FF6 # ⿶
    SURROUND_FROM_LEFT          = 0x2FF7 # ⿷
    SURROUND_FROM_UPPER_LEFT    = 0x2FF8 # ⿸
    SURROUND_FROM_UPPER_RIGHT   = 0x2FF9 # ⿹
    SURROUND_FROM_LOWER_LEFT    = 0x2FFA # ⿺
    OVERLAID                    = 0x2FFB # ⿻

    @staticmethod
    def isIdeoDescription(u):
        if type(u) == str:
            u = ord(u)
        return 0x2FF0 <= u <= 0x2FFB

    def getTransformations(self, upem):
        if self == self.LEFT_TO_RIGHT: # ⿰
            t0 = Identity.scale(0.5, 1.0)
            t1 = Identity.translate(0.5*upem, 0.0).scale(0.5, 1.0)
        elif self == self.ABOVE_TO_BELOW: # ⿱
            t0 = Identity.translate(0.0, 0.5*upem).scale(1.0, 0.5)
            t1 = Identity.scale(1.0, 0.5)
        elif self == self.LEFT_TO_MIDDLE_AND_RIGHT: # ⿲
            t0 = Identity.translate(0.0, 0.0).scale(1/3, 1.0)
            t1 = Identity.translate(upem/3, 0.0).scale(1/3, 1.0)
            t2 = Identity.translate(2*upem/3, 0.0).scale(1/3, 1.0)
            return t0, t1, t2
        elif self == self.ABOVE_TO_MIDDLE_AND_BELOW: # ⿳
            t0 = Identity.translate(0.0, 2*upem/3).scale(1.0, 1/3)
            t1 = Identity.translate(0.0, upem/3).scale(1.0, 1/3)
            t2 = Identity.translate(0.0, 0.0).scale(1.0, 1/3)
            return t0, t1, t2
        elif self == self.FULL_SURROUND: # ⿴
            t0 = Identity
            t1 = Identity.translate(0.1*upem, 0.1*upem).scale(0.8, 0.8)
        elif self == self.SURROUND_FROM_ABOVE: # ⿵
            t0 = Identity
            t1 = Identity.translate(0.1*upem, 0.0).scale(0.8, 0.8)
        elif self == self.SURROUND_FROM_BELOW: # ⿶
            t0 = Identity
            t1 = Identity.translate(0.1*upem, 0.2*upem).scale(0.8, 0.8)
        elif self == self.SURROUND_FROM_LEFT: # ⿷
            t0 = Identity
            t1 = Identity.translate(0.2*upem, 0.1*upem).scale(0.8, 0.8)
        elif self == self.SURROUND_FROM_UPPER_LEFT: # ⿸
            t0 = Identity
            t1 = Identity.translate(0.2*upem, 0.0).scale(0.8, 0.8)
        elif self == self.SURROUND_FROM_UPPER_RIGHT: # ⿹
            t0 = Identity
            t1 = Identity.scale(0.8, 0.8)
        elif self == self.SURROUND_FROM_LOWER_LEFT: # ⿺
            t0 = Identity
            t1 = Identity.translate(0.2*upem, 0.2*upem).scale(0.8, 0.8)
        elif self == self.OVERLAID:
            t0 = Identity
            t1 = Identity
        else:
            assert False
        return t0, t1


# Read database
Hbuild = {}
bases = set()
with open("ids.txt") as f:
    for line in f:
        if count is not None:
            count -= 1
            if not count:
                break
        if line.startswith("#"): continue
        fields = line.split()
        unicode = int(fields[0][2:], 16)
        unichar = fields[1]
        if unicode not in cmap: continue
        build = fields[2].split('[')[0]

        if isCircledNumber(unicode):
            continue

        if unichar == build:
            bases.add(unicode)
            continue

        if any(ord(b) not in cmap and
               not isCircledNumber(b) and
               not IdeoDescription.isIdeoDescription(b)
               for b in build):
            continue

        Hbuild[unicode] = tuple(ord(b) for b in build)

# Prune to ideographs that are not recursively fully available
changed = True
while changed:
    changed = False
    for H,build in list(Hbuild.items()):
        for u in build:
            if not (IdeoDescription.isIdeoDescription(u) or
                    isCircledNumber(u) or
                    u in Hbuild or u in bases):
                print(hex(H), hex(u))
                del Hbuild[H]
                changed = True
                break

print("%d bases; %d ideographs to be matched." % (len(bases), len(Hbuild)))


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
    for u in list(bases) + list(Hbuild) :
        pen = PerContourPen(RecordingPen)
        glyphset[cmap[u]].draw(pen)
        shapes[weight][u] = [recPen.value for recPen in pen.value]


print("Gathering components.")

def recurseBuild(build, t = Identity):
    b = build[0]
    build = build[1:]

    if not IdeoDescription.isIdeoDescription(b):
        if b in bases:
            return [(t, b)], build
        elif isCircledNumber(b):
            return [], build
        else:
            return recurseBuild(Hbuild[b], t)[0], build

    transformations = IdeoDescription(b).getTransformations(upem)

    ret = []
    for trans in transformations:
        trans.transform(t)
        values, build = recurseBuild(build, trans)
        ret.extend(values)

    return ret, build

w0,w1 = WEIGHTS
matches = set()
num_matched = 0
not_matched = 0
mismatch = 0
alternates = defaultdict(list)
HbuildRecursive = {}
Hbuild2 = {}
for H,build in list(Hbuild.items()):

    Hshape0 = shapes[w0][H]
    Hshape1 = shapes[w0][H]

    recursiveBuild,build = recurseBuild(build)
    assert not build
    HbuildRecursive[H] = recursiveBuild

    Rshape = []
    for trans,base in recursiveBuild:
        uShape = transformOutline(trans, shapes[w0][base])
        Rshape.extend(uShape)

    if len(Hshape0) < len(Rshape):
        mismatch += 1
        continue

    matchedOutline0,cost,assignment = matchOutline(Hshape0, Rshape, partial=True)
    if not matchedOutline0:
        not_matched += 1
        continue

    num_matched += 1
    matches.add(H)

    matchedOutline1 = reorderAssignment(shapes[w1][H], assignment)
    assert len(matchedOutline0) == len(matchedOutline1)

    offset = 0;
    rShapes = []
    for trans,base in recursiveBuild:
        baseRshape = shapes[w0][base]
        Rlen = len(baseRshape)
        rShape0 = matchedOutline0[offset:offset+Rlen]
        rShape1 = matchedOutline1[offset:offset+Rlen]
        offset += Rlen

        rShapes.append((base, rShape0, rShape1))
        alternates[base].append((rShape0, rShape1))

    if offset < len(Hshape0):
        otherStrokes0 = matchedOutline0[offset:]
        otherStrokes1 = matchedOutline1[offset:]
        structure0 = outlineStructure(otherStrokes0)

        rShapes.append((structure0, otherStrokes0, otherStrokes1))
        alternates[structure0].append((otherStrokes0, otherStrokes1))

    Hbuild2[H] = rShapes

print("matched: %d not matched: %d mismatch: %d." % (num_matched, not_matched, mismatch))


print("Learning.")
learned = {}
structs = {}
componentDefaultMaster = {}
componentMasters = {}
componentDeltas = {}
componentCoordinates = {}
total_matches = 0
for key,alts in alternates.items():
    total_matches += len(alts)
    if type(key) == int:
        #print("U+%04X: Structure matched %d." % (key, len(alts)))
        pass

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
    #print("Num masters %d max error %d mean-squared error %g" % (k+1, maxError, meanSqError))

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
    #print("Num instances %d num unique instances %d" % (len(instances), len(unique_instances)))
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
    with open("fonts/han/%s/svg/%s.svg" % (serif, "U+%04X" % key if type(key) == int else key), "w") as fd:

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
print("Total matches %d for %d components." % (total_matches, len(alternates)))


print("Building fonts")


style_name = "flat-original-variable"
file_name = "fonts/han/%s/%s-%s.ttf" % (serif,FAMILY_NAME, style_name)
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
file_name = "fonts/han/%s/%s-%s.ttf" % (serif,FAMILY_NAME, style_name)
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
file_name = "fonts/han/%s/%s-%s.ttf" % (serif,FAMILY_NAME, style_name)
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
    if boundsPen.bounds is not None:
        b = [otRound(v) for v in boundsPen.bounds]
    else:
        b = 0, 0, 0, 0
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
