class Ring: #ring of integers

    def __init__(this, num):
        this.num = num

    def is_unit(this):
        return this.num in [-1, 1]

    def is_zero(this):
        return this.num == 0

    def egcd(x, y): #output (g, a1,b1, a2,b2) where a1*x + b1*y = g and a2*x + b2*y = 0, and must satisfy det[[a1,b1],[a2,b2]] = unit
        if y.is_zero(): return (x, Ring.IDENTITY, Ring.ZERO, Ring.ZERO, Ring.IDENTITY)
        g, a1, b1, a2, b2 = Ring.egcd(y, Ring(x.num % y.num))
        a1, b1 = b1, a1-b1*Ring(x.num // y.num)
        a2, b2 = -b2, -a2+b2*Ring(x.num // y.num)
        return (g, a1, b1, a2, b2)

    def format_quotient(this): #format R/(s)
        if this.is_zero():
            return "Z"
        elif this.is_unit():
            return Exception("Why are you formatting the trivial group")
        else:
            return f"Z/{abs(this.num)}Z"

    def __or__(this, that): #HACKY WAY OF SAYING a divides b LMAO
        return (this.is_zero() and that.is_zero()) or (that.num % this.num == 0)

    def __truediv__(this, that):
        return Ring.ZERO if this.is_zero() else Ring(this.num//that.num)

    def __mul__(this, that):
        return Ring(this.num*that.num)

    def __add__(this, that):
        return Ring(this.num+that.num)

    def __sub__(this, that):
        return this+(-that)

    def __neg__(this):
        return Ring(-this.num)

    def __repr__(this):
        return str(this.num)

Ring.ZERO = Ring(0)
Ring.IDENTITY = Ring(1)





"""class Ring: #field of real numbers

    def __init__(this, num):
        this.num = num

    def is_unit(this):
        return this.num != 0

    def is_zero(this):
        return this.num == 0

    def egcd(x, y): #output (g, a1,b1, a2,b2) where a1*x + b1*y = g and a2*x + b2*y = 0, and must satisfy det[[a1,b1],[a2,b2]] = unit
        if x.is_zero(): return (1, Ring.ZERO, Ring.IDENTITY, -Ring.IDENTITY, Ring.ZERO)
        else: return (1, Ring.IDENTITY, Ring.ZERO, -y/x, Ring.IDENTITY)

    def format_quotient(this): #format R/(s)
        if this.is_zero():
            return "R"
        else:
            Exception("Why are you formatting the trivial group")

    def __or__(this, that): #HACKY WAY OF SAYING a divides b LMAO
        return (this.is_zero() and that.is_zero()) or not this.is_zero()

    def __truediv__(this, that):
        return Ring(0) if this.is_zero() else Ring(this.num/that.num)

    def __mul__(this, that):
        return Ring(this.num*that.num)

    def __add__(this, that):
        return Ring(this.num+that.num)

    def __sub__(this, that):
        return this+(-that)

    def __neg__(this):
        return Ring(-this.num)

    def __repr__(this):
        return str(this.num)

Ring.ZERO = Ring(0)
Ring.IDENTITY = Ring(1)"""





def to_ring(M):
    for i,l in enumerate(M):
        M[i] = [Ring(x) for x in l]

def smithform(M): #return matrix, and new basis

    W = len(M[0])
    H = len(M)

    #SMT = resulting M
    S = [[Ring.ZERO]*H for _ in range(H)]
    T = [[Ring.ZERO]*W for _ in range(W)]
    for i in range(H):
        S[i][i] = Ring.IDENTITY
    for i in range(W):
        T[i][i] = Ring.IDENTITY

    #inverse of T
    invT = [[Ring.ZERO]*W for _ in range(W)]
    for i in range(W):
        invT[i][i] = Ring.IDENTITY

    ops = []

    for w,h in zip(range(W,-1,-1), range(H,-1,-1)):
        while not (all(M[i][w-1].is_zero() for i in range(h-1)) and all(M[h-1][j].is_zero() for j in range(w-1))):
            #column ops
            for j in range(w-1):
                x, y = M[h-1][j], M[h-1][j+1]
                g, a1, b1, a2, b2 = Ring.egcd(x, y)
                for i in range(h):
                    M[i][j], M[i][j+1] = a2*M[i][j] + b2*M[i][j+1], a1*M[i][j] + b1*M[i][j+1]
                for i in range(W):
                    T[i][j], T[i][j+1] = a2*T[i][j] + b2*T[i][j+1], a1*T[i][j] + b1*T[i][j+1]
                ops.append((j, j+1, a2, b2, a1, b1))

            #row ops
            for i in range(h-1):
                x, y = M[i][w-1], M[i+1][w-1]
                g, a1, b1, a2, b2 = Ring.egcd(x, y)
                for j in range(w):
                    M[i][j], M[i+1][j] = a2*M[i][j] + b2*M[i+1][j], a1*M[i][j] + b1*M[i+1][j]
                for j in range(H):
                    S[i][j], S[i+1][j] = a2*S[i][j] + b2*S[i+1][j], a1*S[i][j] + b1*S[i+1][j]

            #check if corner number divides everything in its row
            if all(M[h-1][w-1] | M[h-1][j] for j in range(w-1)):
                for j in range(w-1):
                    m = M[h-1][j] / M[h-1][w-1]
                    for i in range(h):
                        M[i][j] -= m*M[i][w-1]
                    for i in range(W):
                        T[i][j] -= m*T[i][w-1]
                    ops.append((j, w-1, Ring.IDENTITY, -m, Ring.ZERO, Ring.IDENTITY))

    for j, j2, a2, b2, a1, b1 in ops[::-1]:
        for i in range(W):
            det = a2*b1 - b2*a1
            invT[i][j], invT[i][j2] = (b1*invT[i][j] - b2*invT[i][j2])/det, (-a1*invT[i][j] + a2*invT[i][j2])/det

    return (M, T, invT)

matrix = [
    [6, -12, 0, -6],
    [0, 6, 0, 3],
    [4, 0, 3, 0],
]

"""matrix = [
    [1,2,3,4],
    [-2,3,-1,5],
    [3,-1,4,-1]
]"""

print("Group presentation of M =")
for l in matrix: print(l)

to_ring(matrix)

matrix, basis, invbasis = smithform(matrix)
S = [Ring.ZERO]*len(matrix[0])
for i,j in zip(range(len(matrix)-1,-1,-1), range(len(matrix[0])-1,-1,-1)): S[j] = matrix[i][j]

factor_list = [s.format_quotient() for s in S if not s.is_unit()]

print()
print("Isomorphism φ : M -> ", " ⊕ ".join(factor_list))

print()
for i,l in enumerate(basis): print(f"φ(e{i+1}) =", [c for c,s in zip(l,S) if not s.is_unit()])

print()
for i,(l,s) in enumerate(zip(invbasis,S)):
    if not s.is_unit(): print(f"φ⁻¹(f{i+1}) =", l)
