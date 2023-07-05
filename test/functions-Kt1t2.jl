# Function to compare two elements of Q(sqrt(2))
function Main.OMFacto.:>(v1::nf_elem, v2::nf_elem)
    v = v1-v2
    cv = representation_matrix(v)[1,:]
    if cv[1]>0 && cv[2]>0 return true end
    if cv[1]<0 && cv[2]<0 return false end
    if cv[1]<0 return cv[1]^2<2*cv[2]^2 end
    return cv[1]^2>2*cv[2]^2
end

function Main.OMFacto.:<(v1::nf_elem,v2::nf_elem) return v2 > v1 end

function Main.OMFacto.:numerator(e::nf_elem) return e*denominator(e) end

import AbstractAlgebra

function Main.OMFacto.:CoeffAndExp(elt::AbstractAlgebra.Generic.LaurentMPolyWrap{fpFieldElem, fpMPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrapRing{fpFieldElem, fpMPolyRing}}, n::nf_elem)
    c=representation_matrix(n)[1,:]
    shift=[0,0]
    # There is no coeff function for Multivariate Laurent Polynomials.
    # We need to shift it via .mpoly but I do not know a way to get the shift done by this function...
    # Therefore, I precompute this shift in a not nice way
    for i in exponent_vectors(elt)
        if i[1] < shift[1] shift[1]=i[1] end
        if i[2] < shift[2] shift[2]=i[2] end
    end
    # Have to make a weird cast here
    ind = [Int64(numerator(c[1]))-shift[1],Int64(numerator(c[2]))-shift[2]]
    return [coeff(elt.mpoly,ind),ind]
end

# Valuation on K[t1,t2] given by v(t1)=1 and v(t2)=sqrt(2)
function Main.OMFacto.:valuation(elt::AbstractAlgebra.Generic.LaurentMPolyWrap{fpFieldElem, fpMPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrapRing{fpFieldElem, fpMPolyRing}})
    # Not efficient at all !!! We are suppose to manipulate positive valuations only ; should be changed if this was not the case anymore
    Rt, z = PolynomialRing(QQ, "z")
    Ka, a = NumberField(z^2-2, "a")
    if elt == 0 return -1 end
    res = -1
    for i in exponent_vectors(elt)
        tmp = i[1]+a*i[2]
        if res == -1
            res = tmp
        elseif res > tmp
            res = tmp
        end
    end
    return res
end

# will not work correctly if we try several K[[t1,...,tn]] cases
function Main.OMFacto.ValueGroup(A::AbstractAlgebra.Generic.LaurentMPolyWrapRing{fpFieldElem, fpMPolyRing})
    z = PolynomialRing(QQ, "z")[2]
    Qa, a = NumberField(z^2-2, "a")
    return [Qa(1),a]
end

function Main.OMFacto.ResidueField(A::AbstractAlgebra.Generic.LaurentMPolyWrapRing{fpFieldElem, fpMPolyRing})
    return base_ring(A) # enough ?
end

# will not work correctly if we try several K[[t1,...,tn]] cases
function Main.OMFacto.GammaCofactors(Gamma::Vector{nf_elem},gamma::nf_elem)
    # Using that I use only two variables here
    M=[representation_matrix(Gamma[1])[1,:];representation_matrix(Gamma[2])[1,:]]
    g=representation_matrix(gamma)[1,:]
    # we want to solve M*v=g
    return M^(-1)*transpose(g) # ca c'est une matrice, je voudrais en faire un vecteur. Neanmoins, res[1] et res[2] ca marche...
end

# In  : F some base field, P an irreducible polynomial above it, s a string to express the new variable.
# Out : The field extension F(P), and an element of this new field that is a root of P.
function Main.OMFacto.IncResField(F,P::fpPolyRingElem,s::String)
    if degree(F)==1 return FiniteField(P,s) end
    p=characteristic(F)
    z=gen(F)
    FF,zz=FiniteField(p,degree(F)*degree(P),s)
    # ICI : I am not sure that elements of F will always convert into elements of FF via FF(z).
    # Should we check that the minimal polynomial of F factorise on FF ?
    y=PolynomialRing(FF,"y")
    tmp=roots(sum([FF(coeff(R,i))*y^i for i in 0:degree(R)])) # not testing if there is a root here...
    return FF, tmp[1]
end
