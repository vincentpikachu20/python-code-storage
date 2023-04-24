# Save audio as 16-bit input.wav
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

def frombytes(B):
    #convert B to an integer
    X = 0
    for b in B[::-1]:
        X = 256*X + b
    return X

def writetowav(S,T):
    #S = list of numbers between -1.0 and 1.0

    with open("input.wav", "rb") as f:
        filebytes = f.read(44 + 16*len(S))

    assert len(S) == len(T), "S and T not same length"

    NumSamples = 2*len(S)
    NumChannels = frombytes(filebytes[22:24])
    assert NumChannels == 2, "input.wav is not stereo"

    SampleRate = 2*frombytes(filebytes[24:28])
    BitsPerSample = frombytes(filebytes[34:36])
    assert BitsPerSample == 16, "input.wav is not 16-bit"

    ByteRate = SampleRate * NumChannels * BitsPerSample//8
    BlockAlign = NumChannels * BitsPerSample//8
    Subchunk2Size = NumSamples * NumChannels * BitsPerSample//8
    ChunkSize = 36 + Subchunk2Size

    fileleftsamples = []
    filerightsamples = []
    for i in range(44,len(filebytes), 4):
        fileleftsamples.append(filebytes[i:i+2])
        filerightsamples.append(filebytes[i+2:i+4])

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

        for i in range(NumSamples):
            #left ear
            if i % 2:
                f.write(fileleftsamples[i//2])
            else:
                x = S[i//2]
                assert -1 <= x <= 1
                X = int(minsize + (maxsize-minsize)*(x+1)/2)
                f.write(tobytes(X, BitsPerSample//8))

            #right ear
            if i % 2:
                f.write(filerightsamples[i//2])
            else:
                x = T[i//2]
                assert -1 <= x <= 1
                X = int(minsize + (maxsize-minsize)*(x+1)/2)
                f.write(tobytes(X, BitsPerSample//8))



# NIKO ascii art

from math import *

asciiart = [
"                                ",
" ##    ##  ##  ##   ##   #####  ",
" ###   ##  ##  ##  ##   ### ### ",
" ####  ##  ##  ## ##    ##   ## ",
" ## ## ##  ##  ####     ##   ## ",
" ##  ####  ##  ## ##    ##   ## ",
" ##   ###  ##  ##  ##   ### ### ",
" ##    ##  ##  ##   ##   #####  ",
"                                "]

S = []
for i in range(len(asciiart[0])):
    if all(a[i]==' ' for a in asciiart):
        S += [0]*6969
    else:
        C = [j for j,a in enumerate(asciiart) if a[i] == '#']
        for j in range(6969):
            r = [min(C), max(C)+1][j % 2]/len(asciiart)
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

# light bulb graph

Ulow = [0.25]*unit
Uhigh = [0.25]*unit

for i in range(unit):
    x = i/unit
    Ulow[i] = -sqrt(1-x**2)
    Uhigh[i] = sqrt(1-x**2)
    if x <= 0.5:
        Ulow[i] = -sqrt(1-x**2) - 0.33

    #squish
    Uhigh[i] = (Uhigh[i]-1)/2.5*1.8+0.8
    Ulow[i] = (Ulow[i]-1)/2.5*1.8+0.8

U = [a for A in zip(Ulow, Uhigh) for a in A]
U = U[::-1] + U

#merge two drawings

T += U
T = [0]*((len(S)-len(T))//2) + T + [0]*((len(S)-len(T)+1)//2) #just to center it

writetowav(S, T)
