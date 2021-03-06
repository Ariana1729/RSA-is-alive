from sage.modules.free_module_integer import IntegerLattice

## Roughly same idea - use lattices to find fac-relations. The issue is that u-vN is actually quite large and has a tiny tiny tiny tiny tiny tiny chance to be smooth haiz
## this was inspired by the old paper

# setup

N = 100000980001501 # given in paper
n = 90
B = N^(1/2)
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
    R[i,i] = round(sqrt(log(p)))
    R[i,n] = round(B*log(p))
rels = []
L = IntegerLattice(R).HKZ()

v = vector([0 for _ in range(n)]+[B*log(N)])
# idk how to cvp in sage so im going to do it manually lol
# now we do lin alg magic
a,b = L.gram_schmidt()
v_cvp = vector(round(i) for i in vector((w*v)/(w*w) for w in a)*b^-1)*L
e_ls = v_cvp[:n]
u,v = 1,1
for e,p in zip(e_ls,p_ls):
    if e>0:
        u*=p^e
    else:
        v*=p^(-e)
# most ideal case
#for i in range(10^6):
#    if i%1000==0:print(i,e_ls)
#    u,v = 1,1
#    for e,p in zip(e_ls,p_ls):
#        e += randint(-3,3)
#        if e>0:
#            u*=p^e
#        else:
#            v*=p^(-e)
#    e = smooth_test(u-v*N)
#    if e != False:
#        break

print(f"N = {N}")
print(f"u = {u}")
print(f"v = {v}")
print(f"u-v*N = {u-v*N}")
print(smooth_test(u-v*N))
