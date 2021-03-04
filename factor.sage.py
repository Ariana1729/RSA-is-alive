

# This file was *autogenerated* from the file factor.sage
from sage.all_cmdline import *   # import sage library

_sage_const_57105727273542262652200264610718642584634617729084089562813 = Integer(57105727273542262652200264610718642584634617729084089562813); _sage_const_24 = Integer(24); _sage_const_1 = Integer(1); _sage_const_2 = Integer(2); _sage_const_0 = Integer(0)
from sage.modules.free_module_integer import IntegerLattice

## setup

N = _sage_const_57105727273542262652200264610718642584634617729084089562813  # 2^200, factors in 3s on my com
n = _sage_const_24 
f = [i+_sage_const_1  for i in range(n)]
h = _sage_const_2 
k = n/_sage_const_2  #idk lol maybe this works?

# generate fac-relations
R = Matrix(ZZ,n+_sage_const_1 ,n+_sage_const_1 )
for i in range(n):
    R[i,n] = round(N*log(Primes()[i]))
R[n,n] = round(N*log(N))
rels = []
for _ in range(n+_sage_const_1 ):
    shuffle(f)
    for i in range(n):
        R[i,i] = N*f[i]
    b = IntegerLattice(R).HKZ()[_sage_const_0 ]
    e = zip([b[i]/R[i,i] for i in range(n)],Primes()[:n])
    u,v = _sage_const_1 ,_sage_const_1 
    for i,j in e:
        if i>_sage_const_0 :
            u*=j**i
        elif i<_sage_const_0 :
            v*=j**(-i)
    print(u,v,u-v*N)
    ### so u-vN isnt smooth what

