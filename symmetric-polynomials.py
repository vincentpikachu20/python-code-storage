#Algorithm currently runs in like O(n!), please optimize

from itertools import permutations

class SymPoly:
    #{(d_1,...,d_n) : coeff}

    def __init__(this, D=None):
        this.poly = {}
        if D != None:
            assert(type(D) == tuple)
            this.poly[tuple(sorted(D)[::-1])] = 1
            this.clear_zeros()

    def clear_zeros(this):
        this.poly = {k:v for k,v in this.poly.items() if v != 0}

    def __add__(A, B):
        C = SymPoly()
        C.poly = dict(A.poly)
        C += B
        return C

    def __iadd__(A, B):
        for k,v in B.poly.items():
            A.poly[k] = A.poly.get(k, 0) + v
        A.clear_zeros()
        return A

    def __neg__(A):
        C = SymPoly()
        C.poly = dict(A.poly)
        for k in C.poly:
            C.poly[k] *= -1
        return C

    def __sub__(A, B):
        return A + (-B)

    def __mul__(A, B):
        C = SymPoly()
        for ka,va in A.poly.items():
            for kb,vb in B.poly.items():
                D = SymPoly._mul_single(ka,kb)
                D = va*vb*D
                C += D
        C.clear_zeros()
        return C

    def __rmul__(A, n):
        C = SymPoly()
        C.poly = dict(A.poly)
        for k in C.poly:
            C.poly[k] *= n
        return C

    def _mul_single(A, B): #extremely inefficient, please optimize
        C = SymPoly()
        Ap = A + (0,)*len(B)
        Bp = B + (0,)*len(A)
        for a in {*permutations(Ap)}:
            for b in {*permutations(Bp)}:
                D = [i+j for i,j in zip(a,b)]
                if sorted(D)[::-1] == D:
                    C += SymPoly(tuple(v for v in D if v != 0))
        return C


class Module: #c1 e1 + c2 e2 + ... just to hash a dict
    def __init__(this, D=None):
        if D == None:
            D = {}
        assert(type(D) == dict)
        this.vector = D

    def __getitem__(this, i):
        return this.vector.get(i)

    def __add__(this, that):
        C = dict(this.vector)
        for k,v in that.vector.items():
            C[k] = C.get(k, 0) + v
        return Module(C)

    def __hash__(this):
        return hash(tuple(sorted(this.vector.items())))

    def __eq__(this, that):
        return hash(this) == hash(that)

    def __repr__(this):
        return " ".join(f"e{k}" + f"^{v}"*(v != 1) for k,v in this.vector.items())


class ElemPoly:
    #{{e_1:d_1, ..., e_k:d_k} : coeff}

    def __init__(this, D=None):
        this.poly = {}
        if D != None:
            assert(type(D) == int)
            if D == 0:
                this.poly[Module()] = 1
            else:
                this.poly[Module({D:1})] = 1

    def clear_zeros(this):
        this.poly = {k:v for k,v in this.poly.items() if v != 0}

    def __add__(A, B):
        C = ElemPoly()
        C.poly = dict(A.poly)
        for k,v in B.poly.items():
            C.poly[k] = C.poly.get(k, 0) + v
        C.clear_zeros()
        return C

    def __neg__(A):
        C = ElemPoly()
        C.poly = dict(A.poly)
        for k in C.poly:
            C.poly[k] *= -1
        return C

    def __sub__(A, B):
        return A + (-B)

    def __mul__(A, B):
        C = ElemPoly()
        for ka,va in A.poly.items():
            for kb,vb in B.poly.items():
                k = ka+kb
                C.poly[k] = C.poly.get(k, 0) + va*vb
        C.clear_zeros()
        return C

    def __rmul__(A, n):
        C = ElemPoly()
        C.poly = dict(A.poly)
        for k in C.poly:
            C.poly[k] *= n
        return C

    def __repr__(this):
        s = ""
        for k,v in this.poly.items():
            if v > 0:
                if s == "": s += f"{v} "*(v != 1) + repr(k)
                else: s += " + " + f"{v} "*(v != 1) + repr(k)
            else:
                if s == "": s += "-" + f"{-v} "*(v != -1) + repr(k)
                else: s += " - " + f"{-v} "*(v != -1) + repr(k)
        return s


mem = {() : ElemPoly(0)}

def findsym(D):
    assert(type(D) == tuple)
    if D in mem: return mem[D]

    D2 = list(D)
    k = 0
    for i in range(len(D)):
        if D[i] == D[0]:
            D2[i] -= 1
            k += 1
        else: break
    D2 = tuple(v for v in D2 if v > 0)
    RHS = SymPoly((1,)*k)*SymPoly(D2) - SymPoly(D)
    ans = ElemPoly(k)*findsym(D2) - sum((v*findsym(k) for k,v in RHS.poly.items()), start=ElemPoly())
    mem[D] = ans
    return ans



# Newton's identities!

for i in range(1,10):
    P = [i]
    print(f"P({list(P)}) = ", findsym(tuple(P)))
