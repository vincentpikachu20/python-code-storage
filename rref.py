# Linear Algebra Done Wrong, chapter 2.2
# reduce to rref and solve the system

from fractions import Fraction

def fractionify(A):
    for i in range(len(A)):
        A[i] = [Fraction(x) for x in A[i]]

def rref(A):
    for i in range(len(A)-1, -1, -1):
        if any(A[i]): # check if nonzero row
            j = A[i].index(1)
            for k in range(i):
                A[k] = [x - A[k][j]*y for x,y in zip(A[k],A[i])]

def echelon(A):
    n = len(A)
    m = len(A[0])
    i = 0 # first 'unlocked' row
    for j in range(m):
        if any(A[k][j] != 0 for k in range(i,n)):
            # find any row containing nonzero at column j
            k = min(k for k in range(i,n) if A[k][j] != 0)
            # swap with row i
            A[i],A[k] = A[k],A[i] # row i is now nonzero
            # reduce
            A[i] = [x/A[i][j] for x in A[i]]
            for k in range(i+1,n):
                A[k] = [x - A[k][j]*y for x,y in zip(A[k],A[i])]
            i += 1

def answer(A):
    n = len(A[0])-1

    print("Matrix form:")
    print(r"$$\begin{pmatrix}")
    for l in A: print(*l[:-1], sep=" & ", end=" \\\\\n")
    print(r"\end{pmatrix} \bf x = \begin{pmatrix}")
    print(*[l[-1] for l in A], sep=r" \\ ")
    print(r"\end{pmatrix}$$")

    print("Vector equation:")
    print("$$")
    print(*["x_{" + str(i+1) + r"} \begin{pmatrix}" + r" \\ ".join(str(l[i]) for l in A) + r"\end{pmatrix}" for i in range(n)], sep=" + ")
    print(r" = \begin{pmatrix}")
    print(*[l[-1] for l in A], sep=r" \\ ")
    print(r"\end{pmatrix}$$")

    fractionify(A)
    echelon(A)
    rref(A)

    #remove zero rows
    while A and not any(A[-1][:-1]):
        if A[-1][-1] != 0:
            print("Solutions: none")
            return
        A.pop()

    pivots = set()

    #find pivots
    for l in A:
        j = min(j for j in range(n) if l[j])
        pivots.add(j)

    free = sorted(set(range(n)) - pivots)

    sol = [0]*n
    var = [{i:0 for i in free} for _ in range(n)]
    for i in free: var[i][i] = 1

    for l in A:
        j = min(j for j in range(n) if l[j])
        sol[j] = l[-1]
        for k in free:
            var[j][k] = -l[k]

    print("Solutions:")
    print(r"$$\bf x = \begin{pmatrix}")
    print(*sol, sep=r" \\ ")
    print(r"\end{pmatrix}")
    for i in free: print(" + x_{" + str(i+1) + r"} \begin{pmatrix}" + r" \\ ".join(str(l[i]) for l in var) + r"\end{pmatrix}")
    if free: print(r",\quad " + ",".join("x_{" + str(i+1) + "}" for i in free) + r"\in\F")
    print("$$")



# exercise 2.2g
A = [
[3,-1,1,-1,2,   5],
[1,-1,-1,-2,-1, 2],
[5,-2,1,-3,3,   10],
[2,-1,0,-2,1,   5]
]

answer(A)
