function reciprocal(F)
    #  In: F some univariate polynomial
    # Out: its reciprocal polynomial
    K=parent(F)
    return K(reverse(F.coeffs))
end

function truncate(F,n)
    #  In: F a univariate polynomial, n a truncation bound
    # Out: F modulo y^n
    if F.length<n+1
        return F
    end
    K=parent(F)
    return K((F.coeffs)[1:n])
end

"""
    AppRoot(F,N)

Computes the N-th approximate root of F

We assume that N divides deg(F)
"""
function AppRoot(F,N)
    #  In: F in K[[x]][y] of degree d, N dividing d
    # Out: the N-th approximate root of F
    d=F.length-1
    if d % N != 0
        error("N must divide the degree of F")
    end
    L=F.parent
    p=Int64(L.base_ring.base_ring.base_ring.n)
    if (d==N) # saving some computations for this easy case
        K=L.base_ring
        return L([invmod(d,p)*F.coeffs[d],K(1)]) # using coeff here gets me a "UndefVarError: coeff not defined"
    end
    m=div(d,N)
    Lz,z=L["zA"] # using PolynomialRing here gets me a "julia undefvarerror: PolynomialRing not defined"
    G=z^N-reciprocal(F)
    psi=L(1)
    inv=L(invmod(N,p))
    k=1
    while k<m+1 # Newton method (diff(G)=N*z^(N-1))
        k*=2
        psi-=truncate(G(psi)*inv,k)
        inv=truncate(2*inv-N*psi^(N-1)*inv^2,k)
    end
    return reciprocal(truncate(psi,m+1))
end

"""
    TaylorExp(F, phi)

Computes the Taylor expansion of F according to phi (see e.g. Modern Computer Algebra, section 9.2)
"""
function TaylorExp(F, phi)
    #  In: F, phi in K[[x]][y]
    # Out: [a_0,...,a_s] s.t. F=a_0+a_1*phi+...+a_s*phi^s
    return _TaylorExp(F,phi)
end

function _TaylorExp(F, phi)
    d=F.length-1
    m=phi.length-1
    k=div(d,m)
    if k<=0
       return [F]
    end
    k=div(k+1,2) # ceil(k/2)
    tmp=phi^k
    q,r=divrem(F,tmp)
    return [_TaylorExp(r,phi);_TaylorExp(q,phi)]
end



"""
    PhiExp(F,Phi)

Computes the Phi-adic expansion of F

Phi is a table of polynomials
"""
function PhiExp(F,Phi)
#  In: F K[[x]][y], Phi a list of such polynomials with "dividing degrees"
# Out: the Phi-adic expansion of F
    k=length(Phi)
    tmp=TaylorExp(F,Phi[k])
    if k>1
        return [PhiExp(i,Phi[1:k-1]) for i in tmp]
    end
    return tmp
end