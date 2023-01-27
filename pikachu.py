def lcm(a,b): #O(log(min(a,b)))
    from math import gcd
    return(a//gcd(a,b)*b)

def factorint(n): #O(sqrt(n))
    a = {}
    while not n % 2:
        if 2 not in a:
            a[2] = 1
        else:
            a[2] = a[2] + 1
        n //= 2
    while not n % 3:
        if 3 not in a:
            a[3] = 1
        else:
            a[3] = a[3] + 1
        n //= 3
    i,r = 5,2
    while i*i <= n:
        if n % i:
            i += r
            r = 6 - r
        else:
            if i not in a:
                a[i] = 1
            else:
                a[i] = a[i] + 1
            n //= i
    if n > 1:
        if n not in a:
            a[n] = 1
        else:
            a[n] = a[n] + 1
    return(a)

def isprime(n): #O(sqrt(n))
    if n <= 1:
        return False
    if n == 2 or n == 3:
        return (True)
    elif not n % 2 or not n % 3:
        return(False)
    i,r = 5,2
    while i * i <= n:
        if n % i == 0:
            return(False)
        i += r
        r = 6 - r
    return(True)

def primesieve(n): #O(n*log(n))
    truth = [True]*(n+1)
    primes = []
    for i in range(2,n+1):
        if truth[i]:
            primes.append(i)
            for j in range(i*i,n+1,i):
                truth[j] = False
    return primes
    
def totient(n): #O(sqrt(n))
    for i in primefactors(n):
        n //= i
        n *= i - 1
    return(n)

def carmichael(n): #O(sqrt(n))
    m = 1
    k = primefactors(n)
    for i in k:
        p = k[i]
        if i == 2 and p >= 3:
            m = 2**(p - 2)
        else:
            m = lcm(m,(i - 1)*i**(p - 1))
    return(m)

# todo: use gaussian elim

def determinant(a,n): #O(n!)
    if n <= 3:
        if n == 3:
            return(a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][0]*a[1][2]*a[2][1] - a[0][1]*a[1][0]*a[2][2] - a[0][2]*a[1][1]*a[2][0])
        elif n == 2:
            return(a[0][0]*a[1][1] - a[0][1]*a[1][0])
        else:
            return(a[0][0])
    else:
        k = 0
        s = 1
        c = []
        for d in a:
            c.append(d[1:])
        for i in range(n):
            m = list(c)
            m.pop(i)
            k += s*a[i][0]*determinant(m,n - 1)
            s *= -1
        return(k)
    
# todo: use gaussian elim
        
def cramer(m,a,n): #O(n*n!)
    d = determinant(m,n)
    if d == 0:
        return([("No solution or infinitely many",1)])
    else:
        k = []
        for i in range(n):
            q = []
            for j in range(n):
                q.append(list(m[j]))
            for j in range(n):
                q[j][i] = a[j]
            k.append(frac(determinant(q,n),d))
        return(k)

def digitcount(n,b=10): #O(log(n))
    c = [0]*b
    while n:
        c[n%b] += 1
        n //= b
    return(c)
    
def digitsum(n,b=10): #O(log(n))
    s = 0
    while n:
        s += n%b
        n //= b
    return(s)

def sqrt2(n): #O(sqrt(n))
    a = b = 1
    i = primefactors(n)
    for p in i:
        a *= p**(i[p]//2)
        b *= 1 + (p - 1)*(i[p]%2)
    return(a,b)

def dictionary(n): #O(1)
    if n == -1:
        a = open("sdcard/qpython/projects3/Pikachu/Words/alphabetical.txt","r")
    else:
        a = open("sdcard/qpython/projects3/Pikachu/Words/bylength.txt","r")
    if n <= 0:
        for i in a.readlines():
            yield(i[:-1])
    else:
        for i in a.readlines():
            if len(i) > n + 1:
                break
            else:
                yield(i[:-1])
    a.close()
    
def diksyunaryo(n): #O(1)
    if n == -1:
        a = open("sdcard/qpython/projects3/Pikachu/Words/alpabeto.txt","r")
    else:
        a = open("sdcard/qpython/projects3/Pikachu/Words/habaan.txt","r")
    if n <= 0:
        for i in a.readlines():
            yield(i[:-1])
    else:
        for i in a.readlines():
            if len(i) > n + 1:
                break
            else:
                yield(i[:-1])
    a.close()

def divisorsum(n): #O(sqrt(n))
    k = 1
    m = primefactors(n)
    for i in m:
        k *= i**(m[i] + 1) - 1
        k //= i - 1
    return(k)

def leastexp(n,m): #O(m)
    from math import gcd
    k = n = n % m
    if gcd(n,m) > 1:
        return(None)
    a = 1
    while k != 1:
        k *= n
        k %= m
        a += 1
    return(a)

def numdivisor(n): #O(sqrt(n))
    m = 1
    p = primefactors(n)
    for i in p:
        m *= p[i] + 1
    return(m)
    
def int2base(x, base): #O(log(n))
    if x < 0:
        sign = -1
    elif x == 0:
        return "0"
    else:
        sign = 1

    x *= sign
    digits = ""
    charmander = "0123456789abcdefghhijklmnopqrstuvwxyz"

    while x:
        digits = charmander[int(x % base)] + digits
        x = x // base

    if sign < 0:
        return("-" + digits)
    else:
        return(digits)

def crt(a,b,c,d): #O(log(b))
    return (c + (a - c)*pow(d,-1,b)*d) % (b*d)

def gcrt(a,b,c,d): #O(log(b))
    g = gcd(b,d)
    if (a-c)%g != 0: return None
    return a%g + crt(a//g,b//g,c//g,d//g)*g
    
def issquare(n): #O(log(n))
    from math import isqrt
    s = isqrt(n)
    return s*s == n

def isresidue(n,p): #O(log(p))
    n %= p
    return n == 0 or pow(n,(p-1)//2,p) == 1

def nonresidue(p):
    k = 2
    while isresidue(k,p):
        k += 1
    return k

def tonelli(n,p):
    S = 0
    Q = p-1
    while Q%2 == 0:
        S += 1
        Q //= 2
    z = nonresidue(p)
    M = S
    c = pow(z,Q,p)
    t = pow(n,Q,p)
    R = pow(n,(Q+1)//2,p)
    while True:
        if t == 0:
            return 0
        elif t == 1:
            return R
        i = 0
        tt = t
        while tt != 1:
            i += 1
            tt = pow(tt,2,p)
        b = pow(c,2**(M-i-1),p)
        M = i
        c = pow(b,2,p)
        t = (t*c) % p
        R = (R*b) % p

def sqrtmod(n,p):
    n %= p
    if n == 0 or p == 2:
        return (n,)
    if not isresidue(n,p):
        return ()
    r = tonelli(n,p)
    r = min(r,p-r)
    return (r,p-r)

def sqrtmodd(a,p,k):
    if k == 1: return sqrtmod(a,p)
    r = []
    for s in sqrtmodd(a,p,k-1):
        if 2*s % p == 0:
            if (a - s**2) % p**k == 0:
                for x in range(0,p**k,p**(k-1)):
                    r.append(x + s)
        else:
            r.append(s +((a-s**2) * pow(2*s,-1,p))%p**k)
    return tuple(sorted(r))

def sqrtmodder(a,n):
    r = [0]
    Q = 1
    for p,e in factorint(n).items():
        P = p**e
        s = sqrtmodd(a,p,e)
        if len(s) == 0: return ()
        nr = []
        for i in s:
            for j in r:
                nr.append(crt(i,P,j,Q))
        r = nr
        Q *= P
    return tuple(sorted(r))

def newbinom(C,n0,r0,n,r,M = None):
    num = 1
    den = 1
    if n0 <= n:
        for i in range(n0+1,n+1):
            num *= i
        for i in range(n0-r0+1,n-r0+1):
            den *= i
    else:
        for i in range(n0,n,-1):
            den *= i
        for i in range(n0-r0,n-r0,-1):
            num *= i
    for i in range(r0+1,r+1):
        den *= i
    for i in range(n-r0,n-r,-1):
        num *= i
    if M == None:
        return (C*num//den)
    else:
        return (C*num*pow(den,-1,M)) % M

def newfact(F,m,n,M = None):
    if m <= n:
        for i in range(m+1,n+1):
            F *= i
        if M == None:
            return F
        else:
            return F % M
    else:
        den = 1
        for i in range(m,n,-1):
            den *= i
        if M == None:
            return F // denom
        else:
            return (F * pow(den,-1,M)) % M

def matmul(a,b,M = None):
    import numpy
    c = numpy.matmul(a,b,dtype=object)
    if M == None:
        return c
    else:
        return c % M

def matpow(a,n,M = None):
    import numpy
    if n == 0:
        return numpy.identity(a.shape[0],dtype=object)
    s = matpow(a,n//2,M)
    ss = matmul(s,s,M)
    if n % 2 == 1:
        return matmul(ss,a,M)
    else:
        return ss

def sequenceplot(f,x0,x1,p = 4):
    from PIL import Image
    terms = [f(x) for x in range(x0,x1+1)]
    y0,y1 = min(terms),max(terms)
    img = Image.new('RGB',(p*(x1-x0+4),p*(y1-y0+4)),(0,0,0))
    pix = img.load()
    for x,y in enumerate(terms):
        for i in range((p+1)//2):
            for j in range((p+1)//2):
                pix[i+p*(x+2),p*(y1-y+2)-j] = (255,255,255)
    img.show()

def icbrt(n):
    l = -1
    r = n+1
    while l+1 < r:
        m = (l+r)//2
        if m**3 == n:
            return m
        if m**3 < n:
            l = m
        else:
            r = m
    return l

def pell(D): #x^2 - D*y^2 = 1
    from math import isqrt
    a0 = isqrt(D)
    m,n = a0,D-a0**2
    if n == 0: #D can't be perfect square
        return None 
    a = (a0+m)//n
    p0,q0 = 1,0
    p,q = a0,1
    while a != 2*a0:
        p0,q0, p,q = p,q, a*p+p0, a*q+q0
        m,n = a*n-m, (D-(a*n-m)**2)//n
        a = (a0+m)//n
    if p**2-D*q**2 == 1: #that's the answer
        return p,q 
    elif p**2-D*q**2 == -1: #square the answer
        return p**2+D*q**2, 2*p*q
    else: print("Wtf is going on in your code")

def genpell(D,N): #x^2 - D*y^2 = n,  when n^2 < D
    from math import isqrt
    a0 = isqrt(D)
    m,n = a0,D-a0**2
    if n == 0: #D can't be perfect square
        return None 
    a = (a0+m)//n
    p0,q0 = 1,0
    p,q = a0,1
    sol = []
    M = p**2-D*q**2
    while M != 1:
        p0,q0, p,q = p,q, a*p+p0, a*q+q0
        m,n = a*n-m, (D-(a*n-m)**2)//n
        a = (a0+m)//n
        M = p**2-D*q**2
        if N % M == 0 and (N > 0) == (M > 0):
            f = isqrt(N//M)
            if f**2 == N//M:
                sol.append((f*p,f*q))
    return sol

def genpell(D,N): #x^2 - D*y^2 = n,  no restriction on n
    t,u = pell(D)
    if N > 0:
        L1 = 0
        L2 = isqrt(N*(t-1)//(2*D))
    else:
        L1 = isqrt((-N)//D)
        L2 = isqrt((-N)*(t+1)//(2*D))
    sol = []
    for y in range(L1,L2+1):
        x2 = N + D*y**2
        x = isqrt(x)
        if x**2 == x2:
            sol.append((x,y))
            sol.append((-x,y))
            if y != 0:
                sol.append((x,-y))
                sol.append((-x,-y))
    return sol

def genfunc(p,q,n,MOD = None): #nth term of the genfunc p(x)/q(x), MAKE SURE ALL COEFFICIENTS WILL BE INTEGERS
    p += [0]*(len(q)-len(p))
    q += [0]*(len(p)-len(q))
    a = []
    for i in range(len(p)):
        s = p[i]
        for j in range(i):
            s -= a[j]*q[i-j]
        a.append(s)
    if n < len(p)-1: return a[n]
    n -= len(p)-1
    M = [[0]*(len(q)-1) for _ in range(len(q)-1)]
    for i in range(len(p)-2):
        M[i+1][i] = 1
    M[0] = [-i//q[0] for i in q[1:]]
    M = numpy.array(M,dtype=object)
    O = numpy.array([[i] for i in a[:0:-1]],dtype=object)
    P = matmul(matpow(M,n,MOD), O, MOD)
    return P[0][0]

def primepi(n):
    r = isqrt(n)
    V = [n//i for i in range(1,r+1)]
    V += list(range(V[-1]-1,0,-1))
    S = {i:i-1 for i in V}
    for p in range(2,i+1):
        if S[p] > S[p-1]:
            for v in V:
                if v < p**2: break
                S[v] -= S[v//p] - S[p-1]
    return S[n]

def primesum(n): #10, 501, 521
    r = isqrt(n)
    V = [n//i for i in range(1,r+1)]
    V += list(range(V[-1]-1,0,-1))
    S = {i:i*(i+1)//2-1 for i in V}
    for p in range(2,r+1):
        print(p)
        if S[p] > S[p-1]:
            for v in V:
                if v < p**2: break
                S[v] -= p*(S[v//p] - S[p-1])
    return S[n]

def totientsieve(n):
    primes = primesieve(n)
    t = list(range(n+1))
    for p in primes:
        for q in range(p,n+1,p):
            t[q] = t[q]//p * (p-1)
    return t

def findprimroot(p):
    if p == 2: return 1
    from sympy import divisors
    d = divisors(p-1)
    for g in range(2,p-1):
        a = g
        good = True
        for i in range(len(d)-2):
            a = (a*pow(g,d[i+1]-d[i],p)) % p
            if a == 1:
                good = False
                break
        if good:
            return g

def floorsum(n,x): # sum of floor(i*x) from i = 1 to n
    if n == 0: return 0
    g = int(x)
    m = int(n*(x-g))
    return g*n*(n+1)//2 + m*n - floorsum(m,1/(x-g))

def floorsum(n,a,b,c): # sum of floor(i*x) from i = 1 to n  where x = (a + sqrt(b))/c
    # CLEAN!!
    if n == 0: return 0
    gc = gcd(a,c)
    if b % gc**2 == 0:
        a //= gc
        b //= gc**2
        c //= gc
    g = (a+isqrt(b))//c
    m = (a*n+isqrt(b*n**2)-g*c*n)//c
    return g*n*(n+1)//2 + m*n - floorsum(m,g*c**2-a*c,b*c**2,b-(a-g*c)**2)





#########
# FENWICK

def getsum(BITTree,i):
    s = 0
    i += 1
    while i > 0:
        s += BITTree[i]
        i -= (i & (-i))
    return s

def updatebit(BITTree, n, i, v):
    i += 1
    while i <= n:
        BITTree[i] += v
        i += i & (-i)

def construct(arr,n):
    BITTree = [0]*(n+1)
    for i in range(n):
        updatebit(BITTree, n, i, arr[i])
    return BITTree


#############




########
# FFT

def fft(A,N,p,glist):
    #want to get [A(1),A(w),A(w^2),...,A(w^(N-1))]
    #N = 1 treeevial
    if N == 1:
        return [A[0] % p]
    #otherwise
    fftAeven = fft(A[::2], N//2, p, glist[::2])*2
    fftAodd = fft(A[1::2], N//2, p, glist[::2])*2
    fftA = [(a+b*g) % p for a,b,g in zip(fftAeven,fftAodd,glist)]
    return fftA

def ifft(fftA,N,p,glist):
    #want to compute the following:
    #F(1) + F(ω) + F(ω^2) + ... + F(ω^(n-1))
    #F(1) + ω^(-1)F(ω) + ω^(-2)F(ω^2) + ... + ω^(1-n)F(ω^(n-1))
    #F(1) + ω^(-2)F(ω) + ω^(-4)F(ω^2) + ... + ω^(2-2n)F(ω^(n-1))
    #...
    #then divide all by N
    A = fft(fftA,N,p,glist)
    A = [A[0]] + A[:0:-1]
    n = pow(N,-1,p)
    A = [a*n % p for a in A]
    return A

def polymul(A,B,n,p = 998244353,g = 3):
    N = 2**n
    #multiply A and B
    #modulo p, (p-1) must be divisible by 2^n
    if (p-1) % N != 0:
        raise ValueError("no illegalll!!!")
    #check if primitive root is ok
    g = pow(g,(p-1)//N,p)
    _g = g
    for _ in range(n):
        if _g == 1:
            raise ValueError("primitive root doesn't make sense")
        _g = (_g*_g) % p
    if _g != 1: 
        raise ValueError("primitive root doesn't make sense")
    #sizes must be <= 2**n
    if len(A) > N or len(B) > N:
        raise ValueError("nooo too big")
    A += [0]*(N - len(A))
    B += [0]*(N - len(B))
    #generate [1,g,g^2,g^3,...]
    glist = [1]
    for _ in range(1,N):
        glist.append(glist[-1]*g % p)
    #time for spicy
    fftA = fft(A,N,p,glist)
    fftB = fft(B,N,p,glist)
    fftC = [a*b % p for a,b in zip(fftA,fftB)]
    C = ifft(fftC,N,p,glist)
    return C

def polymul1e97(A,B,n):
    p = 1000000021000000147887751169
    g = 3
    C = polymul(A,B,n,p,g)
    MOD = 10**9+7
    return [c % MOD for c in C]
