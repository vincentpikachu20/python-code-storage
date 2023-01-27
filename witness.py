from random import randint, choice, shuffle

def randomsol(start,end,w,h): #returns vertices,edges
    finished = False
    while not finished:
        crossed = [[False]*h for _ in range(w)]
        crossed[start[0]][start[1]] = True
        vertices = set(start)
        edges = set()
        curver = start
        while curver != end:
            possible = []
            if curver[0] > 0 and not crossed[curver[0]-1][curver[1]]:
                possible.append((curver[0]-1,curver[1]))
            if curver[0] < w-1 and not crossed[curver[0]+1][curver[1]]:
                possible.append((curver[0]+1,curver[1]))
            if curver[1] > 0 and not crossed[curver[0]][curver[1]-1]:
                possible.append((curver[0],curver[1]-1))
            if curver[1] < h-1 and not crossed[curver[0]][curver[1]+1]:
                possible.append((curver[0],curver[1]+1))
            if len(possible) == 0:
                break
            nextver = choice(possible)
            crossed[nextver[0]][nextver[1]] = True
            edges.add((curver,nextver) if curver < nextver else (nextver,curver))
            vertices.add(nextver)
            curver = nextver
        if curver == end:
            finished = True
    return vertices,edges

def regiondivide(edges,w,h): #returns sectors,color,regionsizes
    sectors = [[-1]*h for _ in range(w)]
    uncolored = True
    color = -1
    regionsizes = []
    while uncolored:
        color += 1
        uncolored = False
        for i in range(w):
            for j in range(h):
                if sectors[i][j] == -1:
                    queue = [(i,j)]
                    uncolored = True
                    break
            if uncolored:
                break
        if not uncolored:
            break
        size = 0
        while len(queue) > 0:
            curcell = queue.pop()
            sectors[curcell[0]][curcell[1]] = color
            size += 1
            if curcell[0] > 0 and sectors[curcell[0]-1][curcell[1]] == -1 and ((curcell[0],curcell[1]),(curcell[0],curcell[1]+1)) not in edges:
                queue.append((curcell[0]-1,curcell[1]))
            if curcell[0] < w-1 and sectors[curcell[0]+1][curcell[1]] == -1 and ((curcell[0]+1,curcell[1]),(curcell[0]+1,curcell[1]+1)) not in edges:
                queue.append((curcell[0]+1,curcell[1]))
            if curcell[1] > 0 and sectors[curcell[0]][curcell[1]-1] == -1 and ((curcell[0],curcell[1]),(curcell[0]+1,curcell[1])) not in edges:
                queue.append((curcell[0],curcell[1]-1))
            if curcell[1] < h-1 and sectors[curcell[0]][curcell[1]+1] == -1 and ((curcell[0],curcell[1]+1),(curcell[0]+1,curcell[1]+1)) not in edges:
                queue.append((curcell[0],curcell[1]+1))
        regionsizes.append(size)
    return [sectors,color,regionsizes]

def tetrispiece(sectors,color,w,h,rotate=False):
    if rotate:
        rotatedir = randint(1,4)
    else:
        rotatedir = 0
    parts = [[(i,j),(w-1-i,h-1-j),(h-i-j,i),(j,w-1-i)][rotatedir] for i in range(w) for j in range(h) if sectors[i][j] == color]
    leftbound = min(i[0] for i in parts)
    upbound = min(i[1] for i in parts)
    cropparts = [(i-leftbound,j-upbound) for i,j in parts]
    dots1 = {(0,0):1,(0,1):2,(0,2):4,(1,0):8,(1,1):16,(1,2):32,(0,3):64,(1,3):128}
    dots2 = {(2,0):1,(2,1):2,(2,2):4,(3,0):8,(3,1):16,(3,2):32,(2,3):64,(3,3):128}
    tet1 = 10240
    tet2 = 10240
    for i,j in cropparts:
        if i <= 1:
            tet1 += dots1[(i,j)]
        else:
            tet2 += dots2[(i,j)]
    return [choice(parts),chr(tet1)+chr(tet2)]

def sunpiece2(sectors,color,w,h):
    parts = [(i,j) for i in range(w) for j in range(h) if sectors[i][j] == color and clues[i][j][0] == 'empty']
    if len(parts) < 2:
        return None,None
    shuffle(parts)
    return parts[0],parts[1]

def negatorpiece(sectors,color,w,h):
    parts = [(i,j) for i in range(w) for j in range(h) if sectors[i][j] == color and clues[i][j][0] == 'empty']
    if len(parts) < 2:
        return None,None
    shuffle(parts)
    return parts[0],parts[1]

w = 7
h = 7

numcolors = 0
edges = []
vertices = []
regionsizes = []

while numcolors < 6 or len(edges) < 20 or 5 not in regionsizes:
    sol = randomsol((0,0),(w-1,h-1),w,h)
    edges = sol[1]
    vertices = sol[0]

    celldivide = regiondivide(edges,w-1,h-1)
    sectors = celldivide[0]
    numcolors = celldivide[1]
    regionsizes = celldivide[2]

cellnums = [(i,j) for i in range(w-1) for j in range(h-1)]
shuffle(cellnums)

clues = [[('empty',)]*h for _ in range(w)]
regionclues = [0]*numcolors

#produce doritos
for _ in range(randint(0,12)):
    i,j = cellnums.pop()
    doritocount = sum([e in edges for e in [((i,j),(i,j+1)),((i+1,j),(i+1,j+1)),((i,j+1),(i+1,j+1)),((i,j),(i+1,j))]])
    if doritocount > 0:
        clues[i][j] = ('dorito',doritocount)
        #regionclues[sectors[i][j]] += 1

#produce colored squares
for _ in range(randint(0,10)):
    i,j = cellnums.pop()
    clues[i][j] = ('color',sectors[i][j])
    regionclues[sectors[i][j]] += 1

#produce hexagons
hexagons = list(vertices)+list(edges)
shuffle(hexagons)
hexagons = hexagons[:randint(0,7)]

#produce broken edges
brokenedges = [((i,j),(i,j+1)) for i in range(w) for j in range(h-1) if ((i,j),(i,j+1)) not in edges]+[((i,j),(i+1,j)) for i in range(w-1) for j in range(h) if ((i,j),(i+1,j)) not in edges]
shuffle(brokenedges)
brokenedges = brokenedges[:randint(0,7)]

#produce tetris blocks
tetriscolor = regionsizes.index(5)
if regionclues[tetriscolor] == 0:
    tetrisblock = tetrispiece(sectors,tetriscolor,w-1,h-1)
    i,j = tetrisblock[0]
    clues[i][j] = ('rtetris',tetrisblock[1])
    #regionclues[sectors[i][j]] += 1

#produce negator
negatorchoices = [i for i in range(0,numcolors-1) if regionsizes[i] >= 2]
if len(negatorchoices) > 0:
    negatorcolor = choice(negatorchoices)
    negatorblock,otherblock = negatorpiece(sectors,negatorcolor,w-1,h-1)
    if negatorblock != None:
        clues[negatorblock[0]][negatorblock[1]] = ('negator',)
        cluetype = choice(['dorito','color','sun'])
        if cluetype == 'dorito':
            clues[otherblock[0]][otherblock[1]] == ('dorito',randint(1,3))
        elif cluetype == 'color':
            clues[otherblock[0]][otherblock[1]] == ('color',randint(0,numcolors-1))
        elif cluetype == 'sun':
            clues[otherblock[0]][otherblock[1]] == ('sun',randint(0,numcolors-1))

#produce suns (paired)
sunchoices = [i for i in range(0,numcolors-1) if regionclues[i] == 0 and regionsizes[i] >= 2]
if len(sunchoices) > 0:
    suncolor = choice(sunchoices)
    sunblock1,sunblock2 = sunpiece2(sectors,suncolor,w-1,h-1)
    if sunblock1 != None:
        clues[sunblock1[0]][sunblock1[1]] = ('sun',suncolor)
        clues[sunblock2[0]][sunblock2[1]] = ('sun',suncolor)

colors = ['\033['+i+';3'+j+'m' for i in '01' for j in '1234567']
shuffle(colors)

showsolution = False

printboard = ''
for j in range(h-1):
    printboard += '⬣' if (0,j) in hexagons else 'o'
    for i in range(w-1):
        if ((i,j),(i+1,j)) in edges:
            printboard += '--' + ('⬣' if ((i,j),(i+1,j)) in hexagons else '-') + '--' + ('⬣' if (i,j) in hexagons else 'o')
        else:
            printboard += '     o' if showsolution else ('     o' if ((i,j),(i+1,j)) in brokenedges else '-----o')
    printboard += '\n'
    sideline = ''
    midline = ''
    for i in range(w-1):
        if ((i,j),(i,j+1)) in edges:
            sideline += '|'
            midline += '⬣' if ((i,j),(i,j+1)) in hexagons else '|'
        else:
            sideline += ' ' if showsolution else (' ' if ((i,j),(i,j+1)) in brokenedges else '|')
            midline += ' ' if showsolution else (' ' if ((i,j),(i,j+1)) in brokenedges else '|')
        if clues[i][j][0] == 'tetris':
            sideline += '     '
            midline += '  ' + clues[i][j][1] + ' '
        elif clues[i][j][0] == 'rtetris':
            sideline += '     '
            midline += ' ↻' + clues[i][j][1] + ' '
        elif clues[i][j][0] == 'stetris':
            sideline += '     '
            midline += ' ¬' + clues[i][j][1] + ' '
        elif clues[i][j][0] == 'rstetris':
            sideline += '     '
            midline += '¬↻' + clues[i][j][1] + ' '
        elif clues[i][j][0] == 'dorito':
            sideline += '     '
            midline += ['     ','  △  ',' △ △ ',' △△△ '][clues[i][j][1]]
        elif clues[i][j][0] == 'negator':
            sideline += '     '
            midline += '  ⅄  '
        elif clues[i][j][0] == 'color':
            sideline += '     '
            midline += '  ' + colors[clues[i][j][1]] + '█\033[0m  '
        elif clues[i][j][0] == 'sun':
            sideline += '     '
            midline += '  ' + colors[clues[i][j][1]] + '☼\033[0m  '
        else:
            sideline += '     '
            midline += '     '
    if ((w-1,j),(w-1,j+1)) in edges:
        sideline += '|'
        midline += '⬣' if ((w-1,j),(w-1,j+1)) in hexagons else '|'
    else:
        sideline += ' ' if showsolution else (' ' if ((w-1,j),(w-1,j+1)) in brokenedges else '|')
        midline += ' ' if showsolution else (' ' if ((w-1,j),(w-1,j+1)) in brokenedges else '|')
    printboard += sideline + '\n' + midline + '\n' + sideline + '\n'
printboard += 'o'
for i in range(w-1):
    if ((i,h-1),(i+1,h-1)) in edges:
        printboard += '--' + ('⬣' if ((i,w-1),(i+1,w-1)) in hexagons else '-') + '--' + ('⬣' if (i,w-1) in hexagons else 'o')
    else:
        printboard += '     o' if showsolution else ('     o' if ((i,h-1),(i+1,h-1)) in brokenedges else '-----o')
print('S'+printboard[1:-1]+'E')
