class Lattice():
    def __init__(self,basis):
        ## We first get a BKZ-reduced basis in sage this does the same as .LLL i believe
        B = Matrix(ZZ,basis).BKZ()
        self.B = B.T # paper uses column sage uses rows haiz
        self.basis = list(B)
        self.dim = len(basis)
        m = len(basis[0]) # ambient space dimension
        self.det = det(B*B.T)
        self.rd = 1 # i have no idea how to compute \gamma_n and i dont think it really matters?
        self.A = sqrt(m)*abs(self.det)^(1/m)# estimate for lambda_1

        ## L_t, pi_t, zeta_t. Note that for my sanity I start t from 0
        self.basis_t = []
        self.det_t = []
        pi_gens_t = None
        self.pi_t = []
        self.zeta_t = []
        for t in range(self.dim):
            self.basis_t.append(list(B[:t]))
            self.det_t.append(det(B[:t]*B[:t].T))
            pi_gens_t = B[:t].right_kernel().matrix().gram_schmidt()[0]
            pi_t = lambda v: sum((w*v)/(w*w)*w for w in pi_gens_t)
            zeta_t = lambda v: v-pi_t(v)
            M = []
            for j in range(m):
                M.append(pi_t(vector(1 if i==j else 0 for i in range(m))))
            self.pi_t.append(Matrix(QQ,m,m,M))
            M = []
            for j in range(m):
                M.append(zeta_t(vector(1 if i==j else 0 for i in range(m))))
            self.zeta_t.append(Matrix(QQ,m,m,M))
    def rho_t(self,t):
        return lambda b:sqrt(self.A^2-(self.pi_t[t]*b).norm()^2)
    def beta_t(self,t):
        return lambda b:self.rho_t(t)(b)^(t-1)*pi^((t-1)/2)/(self.det_t[t]*gamma((t+1)/2))

v = [vector(randint(-99,99) for _ in range(9)) for _ in range(5)]
L = Lattice(v)
pi2 = L.pi_t[2]
zeta2 = L.zeta_t[2]
w = vector(randint(-9,9) for _ in range(9))
print((pi2*w)*(zeta2*w)) # sanity check
print(L.det_t[4]*L.beta_t(4)(w)/L.rho_t(4)(w)^3) # make sure volume of 3-ball is actually 4/3 pi r^3
for i in range(5):
    print(n(L.beta_t(i)(w))) # these should all be between 0 and 1 as they are probabilities

