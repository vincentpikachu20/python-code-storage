from itertools import *
from fractions import Fraction
from copy import deepcopy

class Digit:
    def __init__(this,a):
        this.val = a
        this.repr = str(a)
    def __repr__(this):
        return f'"{this.repr}"'
        
#class Negate:

#class Sqrt:
    
class Concat:
    def __init__(this,*a):
        this.repr = "".join(str(i.val) for i in a)
        this.val = int(this.repr)
    def valid(*a):
        if not all(type(i) == Digit for i in a): return False
        if a[0].val == 0: return False #no leading 0
        return True
    def __repr__(this):
        return f'"{this.repr}"'

class Add:
    def __init__(this,*a):
        this.val = sum(i.val for i in a)
        this.repr = "+".join(f"{i.repr}" for i in a)
    def valid(*a):
        if not all(a[i].val <= a[i+1].val for i in range(len(a)-1)): return False #must be sorted
        if any(type(i) in [Add,Subtract] for i in a): return False #dont do nested addition/subtraction
        return True
    def __repr__(this):
        return f'"{this.repr}"'
        
class Multiply:
    def __init__(this,*a):
        this.val = 1
        for i in a: this.val *= i.val
        this.repr = "*".join(i.repr if type(i) in [Digit,Concat] else f"({i.repr})" for i in a)
    def valid(*a):
        if not all(a[i].val <= a[i+1].val for i in range(len(a)-1)): return False #must be sorted
        if any(type(i) in [Multiply,Divide] for i in a): return False #dont do nested multiplication/division
        return True
    def __repr__(this):
        return f'"{this.repr}"'

class Subtract:
    def __init__(this,a,b):
        this.val = a.val - b.val
        this.repr = f"{a.repr}-{b.repr}"
    def valid(*a):
        if len(a) != 2: return False
        if any(type(i) in [Add,Subtract] for i in a): return False #dont do nested add/subtract
        if a[0].val < a[1].val: return False #only nonnegative
        return True
    def __repr__(this):
        return f'"{this.repr}"'

class Divide:
    def __init__(this,a,b):
        this.val = Fraction(a.val, b.val)
        this.repr = "/".join(i.repr if type(i) in [Digit,Concat] else f"({i.repr})" for i in (a,b))
    def valid(*a):
        if len(a) != 2: return False
        if a[1].val == 0: return False #no division by 0
        if any(type(i) in [Multiply,Divide] for i in a): return False #no nested multiply/dividr
        return True
    def __repr__(this):
        return f'"{this.repr}"'

#class Power:

partitions = []
def genpartitions(S, digs):
    if len(digs) == 0:
        if len(S) == 1: return
        global partitions
        partitions.append(deepcopy(S))
        return
    x = digs.pop()
    for i in range(len(S)):
        S[i].append(x)
        genpartitions(S, digs)
        S[i].pop()
    S.append([x])
    genpartitions(S, digs)
    S.pop()
    digs.append(x)

mem = {}

def enumall(digs, oper1, oper2):
    digs.sort()
    if tuple(digs) in mem: return mem[tuple(digs)]
    ANS = []
    if len(digs) == 1:
        ANS.append(Digit(digs[0]))
    else:
        global partitions
        partitions = []
        genpartitions([], digs)
        party = deepcopy(partitions)
        for P in party:
            for A in product(*[enumall(p, oper1, oper2) for p in P]):
                for S in permutations(A):
                    for Op in oper2:
                        if Op.valid(*S):
                            ANS.append(Op(*S))
    for a in list(ANS):
        for Op in oper1:
            if Op.valid(a):
                ANS.append(a)
    mem[tuple(digs)] = ANS
    return ANS

#INPUT HERE
digs = [1,2,3,4,5]
target = 69

for i in enumall(digs,
   [],
   [Concat, Add, Subtract, Multiply, Divide]):
    if i.val == target: print(target,"=",i.repr)
