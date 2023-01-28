from PIL import ImageDraw, Image
from random import choice

def lockscreengenerate(locksize):
    adjlist = []
    for i in range(locksize**2):
        ipos = []
        iy = i // locksize
        ix = i % locksize
        for j in range(-locksize,locksize):
            for k in [(ix+j,iy),(ix,iy+j),(ix+j,iy+j),(ix+j,iy-j)]:
                if 0 <= k[0] < locksize and 0 <= k[1] < locksize and k != (ix,iy):
                    ipos.append(k[1]*locksize + k[0])
        adjlist.append(ipos)

    available = list(range(locksize**2))
    path = []

    start = choice(available)
    next = start
    opportunity = [1]

    while len(opportunity) > 0:
        if (next - start) % (locksize+1) == 0:
            for i in range(min(start,next),max(start,next)+1,locksize+1):
                available[i] = None
        elif (next - start) % locksize == 0:
            for i in range(min(start,next),max(start,next)+1,locksize):
                available[i] = None
        elif (next - start) % (locksize-1) == 0:
            for i in range(min(start,next),max(start,next)+1,locksize-1):
                available[i] = None
        else:
            for i in range(min(start,next),max(start,next)+1):
                available[i] = None

        path.append(next)
        opportunity = list(set(available).intersection(set(adjlist[next])))
        if len(opportunity) > 0:
            start,next = next,choice(opportunity)

    img = Image.new("RGB",(50*locksize,50*locksize))
    draw = ImageDraw.Draw(img)

    draw.rectangle([(0,0),(50*locksize,50*locksize)],fill=(255,255,255,255))

    for i in range(locksize):
        for j in range(locksize):
            draw.ellipse([(i*50+25-8,j*50+25-8),(i*50+25+8,j*50+25+8)],fill=(66,100,147,255))

    pointlist = []
    for k in path:
        i,j = k//locksize,k%locksize
        pointlist.append((i*50+25,j*50+25))
    draw.line(pointlist,fill=(66,100,147,255),width=5)

    img.save("lockscreen.png")
    img.close()

lockscreengenerate(3)
