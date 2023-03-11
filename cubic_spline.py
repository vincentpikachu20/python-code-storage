from fractions import Fraction as F

#INPUT HERE

t_0 = 1  #starting x term
h = 2    #diff. bet. x terms
n = 5
y = [1, 2, 4, 8, 16, 32]
y = [F(x) for x in y]    #let's use fractions only






from fractions import Fraction as F

class Root3:
    def __init__(this,a,b):
        this.a = a
        this.b = b
    def __mul__(this,that):
        a1,b1 = this.a, this.b
        a2,b2 = that.a, that.b
        return Root3(a1*a2 + 3*b1*b2, a1*b2 + a2*b1)
    def __rmul__(this,that):
        return Root3(that*this.a, that*this.b)   #scalar multiplication
    def __add__(this,that):
        return Root3(this.a + that.a, this.b + that.b)
    def __repr__(this):
        return f"{this.a} + {this.b} √3"

def rootexp(n): # evaluate (-2 + √3)^i for each i = 0...n
    r = [Root3(1,0)] # a + b√3
    for _ in range(n):
        r.append(r[-1] * Root3(-2,1))
    return r



#generate D array, 0 for padding
#write D(x) = Σ 6(y_{n+1} - 2y_n + y_{n-1}) x^n for n = 1,2,...,n-1
D = [0] + [6*(y[i+1] - 2*y[i] + y[i-1]) for i in range(1,n)]

r_i = rootexp(n) # (-2 + √3)^i for each i = 0...n
D_i = sum((D[i] * r_i[i] for i in range(1,n)), start=Root3(0,0)) #evaluate D(-2 + √3)

#solve the system
w_n1 = -D_i.b / r_i[n].b
w_1 = -D_i.a - w_n1 * r_i[n].a

#extract the coeffs of A(x)
#A(x) = x (D(x) + w_1 + w_{n-1} x^n)/(x^2 + 4x + 1)
Num = [F(0),w_1] + D[1:] + [w_n1,F(0)]
w = []

#divide Num(x) by (x^2 + 4x + 1)
for i in range(n+1):
    q = Num[i]
    w.append(q)
    Num[i] -= q
    Num[i+1] -= 4*q
    Num[i+2] -= q

S = []

for i in range(n):
    a = (w[i+1] - w[i])/6
    b = w[i]/2
    c = y[i+1] - y[i] - (w[i+1] + 2*w[i])/6
    d = y[i]
    a,b,c,d = a, b - 3*i*a, c + 3*i**2*a - 2*i*b, d - i**3*a + i**2*b - i*c

    #appropriate scaling for general t and h
    a,b,c,d = a/h**3, b/h**2, c/h, d
    a,b,c,d = a, b - 3*t_0*a, c + 3*t_0**2*a - 2*t_0*b, d - t_0**3*a + t_0**2*b - t_0*c
    S.append((a,b,c,d))

print("Copy-pasteable to Desmos:\n")

print("[" + ",".join(f"({t_0 + i*h},{y[i]})" for i in range(n+1)) + "]")

for i,(a,b,c,d) in enumerate(S):
    print(f"S_{i}(x) = ", end = "")

    #format for printing
    for l,x in [(a,"x^3 "), (b,"x^2 "), (c,"x "), (d," ")]:
        if l != 0:
            print(("+" if l > 0 else "-"), abs(l), x, end="")

    print("\\left\\{", f"{t_0 + i*h} \\le x \\le {t_0 + (i+1)*h}", "\\right\\}")

for i in range(n): print(f"S_{i}'(x)")
for i in range(n): print(f"S_{i}''(x)")






'''
SCRATCH WORK

wlog assume that t_0 = 0 and h = 1

  0     1     2     3
[0,1] [1,2] [2,3] [3,4]

         0  1   2
S =   y0, y1, y2, ..., yn
S' =  z0, z1, z2, ..., zn
S'' = 0,  w1, w2, ..., 0

let S_n(x) -> a (x-n)^3 + b (x-n)^2 + c (x-n) + d

S_n(n) = y_n
S_n'(n) = z_n
S_n''(n) = w_n

d = y_n
a + b + c + d = y_{n+1}

c = z_n
3a + 2b + c = z_{n+1}

2b = w_n
6a + 2b = w_{n+1}

solving for a,b,c,d:

a = (w_{n+1} - w_n)/6
b = w_n/2
c = y_{n+1} - y_n - (w_{n+1} + 2w_n)/6
d = y_n

z_n = y_{n+1} - y_n - (w_{n+1} + 2w_n)/6
z_{n+1} = y_{n+1} - y_n + (2w_{n+1} + w_n)/6

we have a recurrence:
y_{n+1} - y_n - (w_{n+1} + 2w_n)/6 = y_n - y_{n-1} + (2w_n + w_{n-1})/6
6(y_{n+1} - 2y_n + y_{n-1}) = w_{n+1} + 4w_n + w_{n-1}

write D(x) = Σ 6(y_{n+1} - 2y_n + y_{n-1}) x^n for n = 1,2,...,n-1
let A(x) = w_n x^n for n = 0,1,...,n

generating functions:
D(x) = (A(x) - w_1 x)/x + 4A(x) + (A(x) - w_{n-1} x^{n-1})x

A(x) = x (D(x) + w_1 + w_{n-1} x^n)/(x^2 + 4x + 1)

(x^2 + 4x + 1) A(x) = x (D(x) + w_1 + w_{n-1} x^n)

substitute x = -2 ± √3
w_1 + w_{n-1} (-2 - √3)^n = -D(-2 - √3)
w_1 + w_{n-1} (-2 + √3)^n = -D(-2 + √3)

solve for w_1 and w_{n-1}, win!
'''
