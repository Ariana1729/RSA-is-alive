from sage.modules.free_module_integer import IntegerLattice

## Roughly same idea - use lattices to find fac-relations. The issue is that u-vN is actually quite large and has a tiny tiny tiny tiny tiny tiny chance to be smooth haiz
## this was inspired by the old paper

# setup

N = 57105727273542262652200264610718642584634617729084089562813 # 2^200, factors in 3s on my com
n = 30
B = 2^45
p_ls = Primes()[:n]

# smooth test
def smooth_test(x):
    e = []
    for p in p_ls:
        e.append(0)
        while x%p == 0:
            e[-1] += 1
            x = x/p
    return e if x==1 else False

# generate fac-relations
R = Matrix(ZZ,n,n+1)
for i,p in zip(range(n),p_ls):
    R[i,i] = 1
    R[i,n] = round(B*log(p))
rels = []
L = IntegerLattice(R).HKZ()

v = vector([0 for _ in range(n)]+[B*log(N)])
# idk how to cvp in sage so im going to do it manually lol
# now we do lin alg magic
a,b = L.gram_schmidt()
v_cvp = vector(round(i) for i in vector((w*v)/(w*w) for w in a)*b^-1)*L
e_ls = v_cvp[:n]
# most ideal case
for i in range(10^6):
    if i%1000==0:print(i,e_ls)
    u,v = 1,1
    for e,p in zip(e_ls,p_ls):
        e += randint(-3,3)
        if e>0:
            u*=p^e
        else:
            v*=p^(-e)
    e = smooth_test(u-v*N)
    if e != False:
        break

print(u,v,u-v*N)
