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

## IMPORTANT : we assume that a function ConstInv (T, Int64) exists (the user must define it)
# About ConstInv :
# ConstInv(K[[t]],n) is K(1)//n, i.e. base_ring(A)(1)//n
# ConstInv(Qp,n) is either Qp(1//n)=Qp(1)//n or Fp(1//n)=Fp(1)//n ; to be investigated
# Looks like ConstInv always work the same in my examples, so I add it here for the moment ; this might have to be done in "base ring" files later on.
function ConstInv(A,n::Int64)
    return base_ring(A)(1)//n
end


"""
    AppRoot(F::Generic.Poly{T},N::Int) where {T}

Computes the N-th approximate root of F

We assume that N divides deg(F) (an error is raised otherwise)
"""
# This function actually increases the precision on the base field. One should add a parameter for maximum precision ?
# What kind of trunctation ? Powers of pi only ? Valuation truncation ? If valuation truncation, would need a lot more parameters...
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
    inv=ConstInv(A,N)
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
    return _TaylorExp(F,phi,F.length-1)
end

function _TaylorExp(F::Generic.Poly{T}, phi::Generic.Poly{T}, d::Int) where {T}
    # providing the degree is necessary for recursive reminders calls (that can have a smaller degree)
    m=phi.length-1
    k=div(d,m)
    if k<=0
       return [F]
    end
    k=div(k+1,2) # ceil(k/2)
    tmp=phi^k
    q,r=divrem(F,tmp)
    return [_TaylorExp(r,phi,m*k-1);_TaylorExp(q,phi,d-m*k)]
end

"""
    PhiExp(F::Generic.Poly{T},Phi::Vector{Generic.Poly{T}}) where {T}

Computes the Phi-adic expansion of F (i.e. list of lists of ... of lists of elements of the base ring A)

Phi is a table of polynomials ; we assume that Phi[1] has degree 1 (as in NaOl21)
"""
function PhiExp(F::Generic.Poly{T},Phi::Vector{Generic.Poly{T}}) where {T}
#  In: F in A[x], Phi a list of such polynomials with "dividing degrees"
# Out: the Phi-adic expansion of F
    k=length(Phi)
    tmp=TaylorExp(F,Phi[k])
    if k>1
        return [PhiExp(i,Phi[1:k-1]) for i in tmp]
    end
    # Assuming that degree(Phi[1]) is 1 here
    return [coeff(i,0) for i in tmp]
end

# Will be an internal function more than anything
function PhiEval(l::Vector,Phi::Vector{Generic.Poly{T}}) where {T}
#  In: a list from PhiExp(F,phi)
# Out: normally, F
    k=length(Phi)
    L=Phi[1].parent
    if k>1
        res=L(0)
        for i in eachindex(l)
            res+=PhiEval(l[i],Phi[1:k-1])*Phi[k]^(i-1)
        end
        return res
    end
    # k is 1
    res=L(0)
    for i in eachindex(l)
        res+=l[i]*Phi[k]^(i-1)
    end
    return res
end

"""
    PhiExpMonomials(P::Generic.Poly{T},Phi::Vector) where {T}

Computes the Phi-adic expansion of F, expressed as a list of monomials
  (i.e. a list of [c,i_0,...,i_k] providing a coefficient - element of A - and a list of exponents - for the Phi[j]).

Phi is a table of polynomials ; we assume that Phi[1] has degree 1 (as in NaOl21)
"""
function PhiExpMonomials(P::Generic.Poly{T},Phi::Vector) where {T}
    k=length(Phi)
    tmp=TaylorExp(P,Phi[k])
    res=[]
    if k>1
        for i in eachindex(tmp)
            if tmp[i]!=0
                rec=PhiExpMonomials(tmp[i],Phi[1:k-1])
                for j in rec
                    res=[res;[[j;i-1]]]
                end
            end
        end
        return res
    end
    # below we use the assumption that Phi[1] has degree 1 (by taking the 0 coefficient)
    return filter(x->x[1]!=0,[[coeff(tmp[i],0),i-1] for i in eachindex(tmp)])
end

"""
    PhiMonomialsEval(dvt::Vector,Phi::Vector)

Computes the "evaluation" of the Phi-adic expansion (as computed by PhiExpMonomials)
Outputs a polynomial in A[x]
"""
### TODO ? working, but wondering if L should be A and not A[x] when k=1. Does it change anything on an efficiency point of view ?
# In the recursive calls, we are computing several times the same powers of Phi[i] ; i < length(Phi). Could be optimised, but easily ?
function PhiMonomialsEval(dvt::Vector,Phi::Vector)
    L=Phi[1].parent
    if dvt == [] return L(0) end
    k=length(Phi)
    phiexps = unique([i[k+1] for i in dvt])
    d = maximum(phiexps)
    if k>1
        coeffs=[PhiMonomialsEval(filter(x->x[k+1]==i,dvt),Phi[1:k-1]) for i in 0:d]
    else
        tmp=[filter(x->x[k+1]==i,dvt) for i in 0:d]
        coeffs=[i==[] ? L(0) : L(i[1][1]) for i in tmp] # there is a single element in each non empty list.
    end
    # at that points, coeffs is a list of [ai] such that the output will be \sum_i ai*Phi[k]^i
    n=2
    powphi=[Phi[k]] # we will do a divide and conquer approach using powers of 2 (would a binaty decomposition of d be better ? Not sure).
    while n<d
        powphi=[powphi;powphi[end]^2]
        n*=2
    end
    n=1
    while n<=length(powphi)
        tmp = []
        for i in 1:div(length(coeffs),2)
            tmp=[tmp;coeffs[2*i-1]+powphi[n]*coeffs[2*i]]
        end
        if rem(length(coeffs),2) == 1
            tmp=[tmp;coeffs[end]]
        end
        coeffs=tmp
        n+=1
    end
    return L(coeffs[1])
end