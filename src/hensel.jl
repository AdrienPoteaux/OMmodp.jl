function BaseTruncation(elt, n)
    error("BaseTruncation must be defined for our base ring")
end

# In : g and h two coprime polynomial over the current residue field
#      vg and vh the targeted valuations
#      Phi the first key polynomials, vals their valuations
#      Lambda the correcting errors associated to Phi / the base ring
#      e the ramification index (if we had a single slope, g and h are given in the "classical" form
#        - we will have to make g(y^e) in the end, we save some computations for gcdx ; if there were >1 slopes, e is 1)
# Out : representants for g, h and their cofactors.
## g and h are not Generic.Poly (on examples they are directly handled by Flint and not AbstractAlgebra)
function PhiInitHensel(g,h,vg,vh,Phi::Vector,vals::Vector,Lambda::Vector,e::Int)
    y=gen(parent(g))
    s,t=gcdx(g,h)[2:3]
    if (e>1)
        g=g(y^e)
        h=h(y^e)
        s=s(y^e)
        t=t(y^e)
    end
    G=Representant(g,vg,Phi,vals,Lambda)
    H=Representant(h,vh,Phi,vals,Lambda)
# ATTENTION : S et T ont des valuations negatives... donc n'appartiennent pas forcément au monde de base
    S=Representant(s,-vg,Phi,vals,Lambda)
    T=Representant(t,-vh,Phi,vals,Lambda)
    return G,H,S,T
end

# Return the truncation of dvt (any element of valuation greater than n is removed)
### WARNING : dvt is probably a pointer, so we are actually working in place here
function _PhiTruncate(dvt,vals::Vector,n)
    res=dvt
    k=length(vals)
    if k==0
        return BaseTruncation(dvt,n) # user must define it for our base ring.
    end
    for i in eachindex(dvt)
        res[i]=_PhiTruncate(dvt[i],vals[1:k-1],n-(i-1)*vals[k])
    end
    return res
end

function PhiTruncate(P::Generic.Poly{T}, Phi::Vector, vals::Vector, n) where {T}
    dvt=PhiExp(P,Phi)
    dvt=_PhiTruncate(dvt,vals,n)
    return PhiEval(dvt,Phi)
end

# In  : mu(F-G*H)>=n+mu(F), mu(S*G+T*H-1)>=n
# Out : n becomes 2*N
# vG and vH provide the valuations of respectively G and H
# ATTENTION : S et T ont des valuations negatives... donc n'appartiennent pas forcément au monde de base
function PhiHenselStep(F::Generic.Poly{X},G::Generic.Poly{X},H::Generic.Poly{X},S::Generic.Poly{X},T::Generic.Poly{X},Phi::Vector,vals::Vector,n,vG,vH) where {X}
    E=PhiTruncate(F-G*H,Phi,vals,2*n+vG+vH)
    q,r=divrem(S*E,H)
    q=PhiTruncate(q,Phi,vals,vG+vH+2*n)
    r=PhiTruncate(r,Phi,vals,vG+vH+2*n)
    Gt=PhiTruncate(G+E*T+q*G,Phi,vals,vG+2*n)
    Ht=PhiTruncate(H+r,Phi,vals,vH+2*n)
    B=PhiTruncate(S*Gt+T*Ht-1,Phi,vals,2*n)
    c,d=divrem(S*B,Ht) ## integer division error ici ?
    c=PhiTruncate(c,Phi,vals,2*n)
    d=PhiTruncate(d,Phi,vals,2*n)
    St=PhiTruncate(S-d,Phi,vals,2*n-vG)
    Tt=PhiTruncate(T-B*T-c*Gt,Phi,vals,2*n-vH)
    return Gt,Ht,St,Tt
end

# In  : F in A[x], g and h two residual polynomials (modified if e = 1), vg and vh the "target" valuations (mu(F)=vg+vh in particular)
#       Phi, vals and Lambda resp. the lists of key polynomials, valuations and correcting terms.
#       e the last partial ramification index (or 1 if we are using modified residual polynomials)
#       n the precision we want for the factorisation
# Out : G and H such that mu(F-GH)>mu(F)+n
### Should we also output S and T ?
## g and h are not Generic.Poly (on examples they are directly handled by Flint and not AbstractAlgebra)
function PhiHensel(F::Generic.Poly{Y},g,h,vg,vh,Phi::Vector,vals::Vector,Lambda::Vector,e::Int,n) where {Y}
    G,H,S,T=PhiInitHensel(g,h,vg,vh,Phi,vals,Lambda,e)
    k=PhiVal(PhiExp(F-G*H,Phi),vals)-vg-vh
    tmp=PhiVal(PhiExp(S*G+T*H-1,Phi),vals) # not a direct minimum so that overloading < for nfelem is sufficient (and not <= too)
    if tmp<k k=tmp end # not sure if this minimal is required or if we cannot directly know that from vals ?
    while k < n
        G,H,S,T=PhiHenselStep(F,G,H,S,T,Phi,vals,k,vg,vh)
        k*=2
    end
    return G,H
end