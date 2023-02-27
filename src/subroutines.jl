function reciprocal(F::Generic.Poly{T}) where{T}
    #  In: F some univariate polynomial
    # Out: its reciprocal polynomial
    K=parent(F)
    return K(reverse(F.coeffs))
end

function truncate(F::Generic.Poly{T},n::Int) where {T}
    #  In: F a univariate polynomial, n a truncation bound
    # Out: F modulo y^n
    if F.length<n+1
        return F
    end
    K=parent(F)
    return K((F.coeffs)[1:n])
end

"""
    AppRoot(F::Generic.Poly{T},N::Int) where {T}

Computes the N-th approximate root of F

We assume that N divides deg(F)

IMPORTANT : we assume that a function ConstInv (T, Int64) exists (the user must define it)
"""
# About ConstInv :
# ConstInv(K[[t]],n) is K(1)//n, i.e. A.base_ring(1)//n
# ConstInv(Qp,n) is either Qp(1//n)=Qp(1)//n or Fp(1//n)=Fp(1)//n ; to be investigated
function AppRoot(F::Generic.Poly{T},N::Int) where {T}
    #  In: F in A[x] of degree d, N dividing d
    # Out: the N-th approximate root of F
    d=F.length-1
    if d % N != 0
        error("N must divide the degree of F")
    end
    L=parent(F)
    A=base_ring(L)
    if (d==N) # saving some computations for this easy case
        return L([ConstInv(A,d)*F.coeffs[d],A(1)]) # using coeff here gets me a "UndefVarError: coeff not defined"
    end
    m=div(d,N)
    Lz,z=L["zA"] # using PolynomialRing here gets me a "julia undefvarerror: PolynomialRing not defined"
    G=z^N-reciprocal(F)
    psi=L(1)
    inv=ConstInv(A,N) ## La aussi, je pense qu'un 1//K(N) doit le faire en mieux SOUCI GENERIQUE
    k=1
    while k<m+1 # Newton method (diff(G)=N*z^(N-1))
        k*=2
        psi-=truncate(G(psi)*inv,k)
        if k<m+1
            inv=truncate(2*inv-N*psi^(N-1)*inv^2,k)
        end
    end
    return reciprocal(truncate(psi,m+1))
end

"""
    TaylorExp(F::Generic.Poly{T}, phi::Generic.Poly{T}) where {T}

Computes the Taylor expansion of F according to phi (see e.g. Modern Computer Algebra, section 9.2)
"""
function TaylorExp(F::Generic.Poly{T}, phi::Generic.Poly{T}) where {T}
    #  In: F, phi in A[x]
    # Out: [a_0,...,a_s] s.t. F=a_0+a_1*phi+...+a_s*phi^s
    return _TaylorExp(F,phi)
end

function _TaylorExp(F::Generic.Poly{T}, phi::Generic.Poly{T}) where {T}
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
    PhiExp(F::Generic.Poly{T},Phi::Vector{Generic.Poly{T}}) where {T}

Computes the Phi-adic expansion of F

Phi is a table of polynomials
"""
function PhiExp(F::Generic.Poly{T},Phi::Vector{Generic.Poly{T}}) where {T}
#  In: F in A[x], Phi a list of such polynomials with "dividing degrees"
# Out: the Phi-adic expansion of F
    k=length(Phi)
    tmp=TaylorExp(F,Phi[k])
    if k>1
        return [PhiExp(i,Phi[1:k-1]) for i in tmp]
    end
    return tmp
end