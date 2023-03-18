from fractions import Fraction
from math import factorial

#TO-DO: implement FFT

class GenFunc:
    def __init__(this, A, maxdeg = 100):
        this.poly = list(A)
        this.maxdeg = maxdeg
        this.truncate()

    def copy(this):
        return GenFunc(this.poly, this.maxdeg)

    def truncate(this, maxdeg = None):
        if maxdeg != None:
            this.maxdeg = min(this.maxdeg, maxdeg)
        for _ in range(len(this.poly) - this.maxdeg-1):
            this.poly.pop()
        this.reducezeros()

    def reducezeros(this):
        while len(this.poly) > 0 and this.poly[-1] == 0:
            this.poly.pop()

    def frac(a, b):
        if type(n) == int: return n
        if n.denominator == 1: return n.numerator
        return n

    def __repr__(this):
        return f"GenFunc([{', '.join(str(i) for i in this.poly)}], maxdeg = {this.maxdeg})"

    def __getitem__(this, i):
        if i >= len(this.poly): return 0
        return this.poly[i]

    def __setitem__(this, i, v):
        if i > this.maxdeg: return
        if v != 0:
            while len(this.poly) <= i:  #if v == 0, no need to extend list
                this.poly.append(0)
            this.poly[i] = v
        elif i < len(this.poly):
            this.poly[i] = 0
            this.reducezeros()

    def __eq__(this, that):
        for i in range(min(this.maxdeg, that.maxdeg)):
            if this[i] != that[i]:
                return False
        return True

    def __iadd__(this, that):
        this.truncate(that.maxdeg)
        for i in range(this.maxdeg+1):
            if i >= len(that.poly): break
            this[i] += that[i]
        this.reducezeros()
        return this

    def __add__(this, that):
        A = this.copy()
        A += that
        return A

    def __isub__(this, that):
        this.truncate(that.maxdeg)
        for i in range(this.maxdeg+1):
            if i >= len(that.poly): break
            this[i] -= that[i]
        this.reducezeros()
        return this

    def __sub__(this, that):
        A = this.copy()
        A -= that
        return A

    def __neg__(this):
        A = this.copy()
        for i in range(len(A.poly)):
            A[i] = -A[i]
        return A

    def __imod__(this, m):
        for i in range(len(this.poly)):
            this[i] %= m
        this.reducezeros()
        return this

    def diff(this, n=1):
        assert n >= 0
        if n == 0:
            return this.copy()

        A = GenFunc([], this.maxdeg - n)
        f = factorial(n)
        for i in range(n, len(this.poly)):
            A.poly.append(this[i]*f)
            f //= i-n+1
            f *= i+1
        return A

    def xdiff(this, n=1):
        assert n >= 0
        if n == 0:
            return this.copy()

        A = GenFunc([], this.maxdeg)
        for i in range(len(this.poly)):
            A.poly.append(this[i]*i**n)
        return A

    def integrate(this, n=1):
        assert n >= 0
        if n == 0:
            return this.copy()

        A = GenFunc([], this.maxdeg + n)
        A.poly = [0]*n
        f = factorial(n)
        for i in range(len(this.poly)):
            A.poly.append(Fraction(this[i], f))
            f //= i+1
            f *= i+n+1
        return A

    def xintegrate(this, n=1):
        assert n >= 0
        if n == 0:
            return this.copy()

        A = GenFunc([], this.maxdeg)
        for i in range(len(this.poly)):
            A.poly.append(Fraction(this[i], (i+1)**n))
        return A

    def __rmul__(this, n):
        if n == 0:
            return GenFunc([], this.maxdeg)

        A = this.copy()
        for i in range(len(this.poly)):
            A[i] *= n
        return A

    def __mul__(this, that): #runs in O(n^2). TO-DO: implement FFT
        if type(that) != GenFunc:
            return that * this

        A = GenFunc([], min(this.maxdeg, that.maxdeg))
        for i in range(len(this.poly)):
            for j in range(len(that.poly)):
                if i+j > A.maxdeg: break
                A[i+j] += this[i]*that[j]
        A.reducezeros()
        return A

    def __rtruediv__(this, that):
        return GenFunc([that]) / this

    def __truediv__(this, that): #runs in O(n^2). divide using long division. is there a faster algorithm?
        if type(that) != GenFunc:
            return this * Fraction(1, that)

        if len(that.poly) == 0:
            raise ZeroDivisionError("division by zero")

        # check if it'll give a valid genfunc
        k = 0
        while this[k] == 0 and that[k] == 0:
            k += 1
        if that[k] == 0:
            raise ZeroDivisionError("division result is not a valid power series")

        Q = GenFunc([], min(this.maxdeg, that.maxdeg) - k)
        D = this.copy()
        D.truncate(that.maxdeg)
        for i in range(k, D.maxdeg+1):
            x = Fraction(D[i], that[k])
            Q.poly.append(x)

            allzeros = True
            for j in range(k, len(that.poly)):
                if i+j-k > D.maxdeg: break
                D[i+j-k] -= that[j]*x
                allzeros &= (D[i+j-k] == 0)

            if allzeros: break # terminate if it's all zeros
        Q.reducezeros()
        return Q

    def __pow__(this, n): #binary exp is fast enough
        if n < 0:
            return 1/this**(-n)
        if n == 0:
            return GenFunc([1], this.maxdeg)
        if n == 1:
            return this.copy()
        A = this**(n//2)
        A *= A
        if n % 2 == 1:
            A *= this
        return A

    def __call__(this, that): #composition of generating functions
        if type(that) != GenFunc:
            raise TypeError("should only compose with GenFunc")
        A = GenFunc([], min(this.maxdeg, that.maxdeg))
        X = GenFunc([1], A.maxdeg)
        for i in range(len(this.poly)):
            if i > A.maxdeg:
                break
            A += X*this[i]
            X *= that
        return A

    def inverse(this): #g such that f(g(x)) = x, this is super slow, O(n^3), is there a better way?
        if this[0] != 0 or this[1] == 0:
            raise Exception("inverse only exists iff f(0) = 0 and f'(0) != 0")

        G = GenFunc([0, Fraction(1, this[1])], 1)
        for i in range(2, this.maxdeg+1):
            G.maxdeg += 1

            #r = [-, g, g^2, g^3, ...]
            R = [0, G]
            for _ in range(2, i+1):
                R.append(R[-1] * G)

            c = sum(this[j]*R[j][i] for j in range(2, i+1))
            G[i] = -Fraction(c, this[1])

        G.reducezeros()
        return G

# These next two functions are a bit weird to express

def log(F):
    if F[0] == 0:
        raise ZeroDivisionError("log 0 is undefined")
    G = (F.diff()/F).integrate()
    if F[0] == 1:
        return G
    else:
        G[0] = f"log {F[0]}"   # the answer would be a bit weird
        return G

def exp(F):
    G = GenFunc([1], F.maxdeg)
    H = GenFunc([1], F.maxdeg)
    dF = F.diff()
    f = 1
    for i in range(1, F.maxdeg+1):
        H = H.diff() + H*dF
        f *= i
        G[i] = Fraction(H[0], f)
    if F[0] == 0:
        return G
    else:
        from sympy import symbols
        e = symbols('e')  # yeah it needs e at this point... but only if e^f != 1
        return e**F[0] * G
    

# Examples

A = GenFunc([0,1,2,3], 15)
B = GenFunc([4,5,6,7])

print(A + B)
print(A * B)
print(A - B)
print(A / B)
print(A ** 3)
print(B ** -1)
print(A(B))
print(A.inverse())
print(A.diff())
print(A.integrate())
