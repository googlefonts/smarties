from fontTools.fontBuilder import FontBuilder
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.pens.cu2quPen import Cu2QuPen, Cu2QuMultiPen
from fontTools.ttLib.tables._g_l_y_f import GlyphCoordinates
from fontTools.ttLib.tables.TupleVariation import TupleVariation

def createFontBuilder(font, family_name, style, chars, extraGlyphs=[]):
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
        familyName=dict(en=family_name),
        styleName=dict(en=style),
    )

    fb = FontBuilder(upem, isTTF=True)
    fb.setupHead(unitsPerEm=upem, created=font['head'].created, modified=font['head'].modified)
    fb.setupNameTable(nameStrings)
    fb.setupGlyphOrder(subset_glyphOrder)
    fb.setupCharacterMap(subset_cmap)
    fb.setupHorizontalMetrics(subset_metrics)
    hhea = font['hhea']
    fb.setupHorizontalHeader(ascent=hhea.ascent, descent=hhea.descent)
    os2 = font['OS/2']
    fb.setupOS2(sTypoAscender=os2.sTypoAscender, usWinAscent=os2.usWinAscent, usWinDescent=os2.usWinDescent)
    fb.setupPost(keepGlyphNames=False)

    return fb

def createTTGlyphPen():
    return TTGlyphPen(None, outputImpliedClosingLine=True)
def createCu2QuPen(pen):
    # reverse_direction=True was broken in
    # https://github.com/fonttools/fonttools/pull/2995
    return Cu2QuPen(pen, 1, reverse_direction=False)
def createCu2QuMultiPen(pens):
    # reverse_direction=True was broken in
    # https://github.com/fonttools/fonttools/pull/2995
    return Cu2QuMultiPen(pens, 1, reverse_direction=False)

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

def setupVariableFont(glyphSets, weights):
    assert len(weights) == 2
    variations = {}
    tag = 'wght'
    glyphs = glyphSets[weights[0]]
    varGlyphs = glyphSets[weights[1]]
    for glyphName in glyphs:
        glyph = glyphs[glyphName]
        varGlyph = varGlyphs[glyphName]

        coords = getCoords(varGlyph) - getCoords(glyph)

        coords.extend([(0,0), (0,0), (0,0), (0,0)]) # Phantom points TODO
        axes = {tag: (0, 1, 1)}
        tv = TupleVariation(axes, coords)
        variations[glyphName] = [tv]

    axes = [(tag, weights[0], weights[0], weights[1], tag)]
    return glyphs, variations, axes

def fixLsb(fb):
    metrics = fb.font["hmtx"].metrics
    glyf = fb.font["glyf"]
    for glyphname in glyf.keys():
        v = getattr(glyf.glyphs[glyphname], "xMin", 0)
        metrics[glyphname] = (metrics[glyphname][0], v)
