from sage.modules.free_module_integer import IntegerLattice

## setup

N = 57105727273542262652200264610718642584634617729084089562813 # 2^200, factors in 3s on my com
n = 24
f = [i+1 for i in range(n)]
h = 2
k = n/2 #idk lol maybe this works?

# generate fac-relations
R = Matrix(ZZ,n+1,n+1)
for i in range(n):
    R[i,n] = round(N*log(Primes()[i]))
R[n,n] = round(N*log(N))
rels = []
for _ in range(n+1):
    shuffle(f)
    for i in range(n):
        R[i,i] = N*f[i]
    b = IntegerLattice(R).HKZ()[0]
    e = zip([b[i]/R[i,i] for i in range(n)],Primes()[:n])
    u,v = 1,1
    for i,j in e:
        if i>0:
            u*=j^i
        elif i<0:
            v*=j^(-i)
    print(u,v,u-v*N)
    ### so u-vN isnt smooth ...

