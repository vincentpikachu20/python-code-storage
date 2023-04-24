# Open the resulting output.wav in Audacity



def tobytes(X, n):
    #write X as n-byte string
    #.wav writes integers backwards for some reason

    #tobytes(0x1234, 2) -> 0x3412
    B = []
    for _ in range(n):
        B.append(X % 256)
        X //= 256
    return bytes(B)

def writetowav(S):
    #S = list of numbers between -1.0 and 1.0

    NumSamples = len(S)
    NumChannels = 1
    SampleRate = 44100
    BitsPerSample = 16

    ByteRate = SampleRate * NumChannels * BitsPerSample//8
    BlockAlign = NumChannels * BitsPerSample//8
    Subchunk2Size = NumSamples * NumChannels * BitsPerSample//8
    ChunkSize = 36 + Subchunk2Size

    with open("output.wav", "wb") as f:

        f.write(b"RIFF") #ChunkID
        f.write(tobytes(ChunkSize, 4)) #ChunkSize
        f.write(b"WAVE") #Format

        f.write(b"fmt ") #Subchunk1ID
        f.write(tobytes(16, 4)) #Subchunk1Size
        f.write(tobytes(1, 2)) #AudioFormat

        f.write(tobytes(NumChannels, 2)) #NumChannels
        f.write(tobytes(SampleRate, 4)) #SampleRate
        f.write(tobytes(ByteRate, 4)) #ByteRate
        f.write(tobytes(BlockAlign, 2)) #BlockAlign
        f.write(tobytes(BitsPerSample, 2)) #BitsPerSample

        f.write(b"data") #Subchunk2ID
        f.write(tobytes(Subchunk2Size, 4))

        minsize = -2**(BitsPerSample-1)
        maxsize = 2**(BitsPerSample-1) - 1
        for x in S:
            assert -1 <= x <= 1
            X = int(minsize + (maxsize-minsize)*(x+1)/2)
            f.write(tobytes(X, BitsPerSample//8))



# NIKO ascii art

from math import *
from random import *

asciiart = [
"                               ",
" ##    ##  ##  ##   ##  #####  ",
" ###   ##  ##  ##  ##  ### ### ",
" ####  ##  ##  ## ##   ##   ## ",
" ## ## ##  ##  ####    ##   ## ",
" ##  ####  ##  ## ##   ##   ## ",
" ##   ###  ##  ##  ##  ### ### ",
" ##    ##  ##  ##   ##  #####  ",
"                               "]

S = []
for i in range(len(asciiart[0])):
    if all(a[i]==' ' for a in asciiart):
        S += [0]*6969
    else:
        C = [j for j,a in enumerate(asciiart) if a[i] == '#']
        for j in range(6969):
            r = choice(C)
            r = (r + uniform(0,1))/(len(asciiart))
            S.append(1 - 2*r)


# Niko graph

unit = 6969*2
Tlow = [0.25]*unit*2
Thigh = [0.25]*unit*2

# 1-0.3x^2 {0 <= x <= 1.6}
# 0.1x^2 {0 <= x <= 1.6}
# 9(x-0.69)^2-1 {0.69<=x<=1}
# sqrt(x)+0.5 {0<=x<=0.9}
# 1.2-0.1x^2 {0<=x<=0.5}
# -sqrt(1-x^2) {0<=x<=1}

for i in range(unit*2):
    x = i/unit
    if 0 <= x <= 1.6:
        if 1-0.3*x**2 >= 0.1*x**2:
            Thigh[i] = 1-0.3*x**2
            Tlow[i] = 0.1*x**2
    if 0.69 <= x <= 1:
        Tlow[i] = 9*(x-0.69)**2-1
    if 0 <= x <= 0.9:
        Thigh[i] = max(Thigh[i], sqrt(x)+0.5)
    if 0 <= x <= 0.5:
        Thigh[i] = max(Thigh[i], 1.2-0.1*x**2)
    if 0 <= x <= 1:
        Tlow[i] = min(Tlow[i], -sqrt(1-x**2))

    #squish
    Thigh[i] = (Thigh[i]+1)/2.5*1.8-0.9
    Tlow[i] = (Tlow[i]+1)/2.5*1.8-0.9

T = [a for A in zip(Tlow, Thigh) for a in A]
T = T[::-1] + T

writetowav(S + T)
