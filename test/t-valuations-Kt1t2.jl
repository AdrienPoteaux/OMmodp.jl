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

# exact same function over Fq -> mutualiser avec du Union dans les titres ?
function Main.OMFacto.:CoeffAndExp(elt::AbstractAlgebra.Generic.LaurentMPolyWrap{fqPolyRepFieldElem, fqPolyRepMPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrapRing{fqPolyRepFieldElem, fqPolyRepMPolyRing}}, n::nf_elem)
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

# exact same function over Fq
function Main.OMFacto.:valuation(elt::AbstractAlgebra.Generic.LaurentMPolyWrap{fqPolyRepFieldElem, fqPolyRepMPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrapRing{fqPolyRepFieldElem, fqPolyRepMPolyRing}})
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

# same over Fq
function Main.OMFacto.ValueGroup(A::AbstractAlgebra.Generic.LaurentMPolyWrapRing{fqPolyRepFieldElem, fqPolyRepMPolyRing})
    z = PolynomialRing(QQ, "z")[2]
    Qa, a = NumberField(z^2-2, "a")
    return [Qa(1),a]
end

# same over Fq
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

# la meme sur Fq
function Main.OMFacto.IncResField(F,P::fqPolyRepPolyRingElem,s::String)
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

function TestPhiValKt1t2()
    R, (t1, t2) = LaurentPolynomialRing(GF(149), ["t1", "t2"])
    A, x = PolynomialRing(R, "x")
    phi1=x^12-A(5*t1^3*t2^4)#v(x) = 1/4+sqrt(2)/3
    phi2=phi1^4-2^9*3^8*t1^11*t2^7*x^18# v(phi1) = 31/8+13/4*sqrt(2)
    F=phi2^3-A(144*3^25*5^14*t1^38*t2^24)# v(phi2) =  38/3+8*sqrt(2) ; 144 = 2^(-858) mod 149
# This is a big irreducible example I used in our exchange about multivariate Puiseux expansions.
# I assume the reduction modulo 149 is fine.

    Rt, z = PolynomialRing(QQ, "z")
    Ka, a = NumberField(z^2-2, "a")

    Phi = [x,phi1,phi2]
    vals = [1//4+a//3,31//8+13//4*a,38//3+8*a]

    elt = PhiExp(F,Phi)
    return PhiVal(elt,vals) == 38+24*a
end

function TestPhiNewtonPolygonKt1t2()
    R, (t1, t2) = LaurentPolynomialRing(GF(149), ["t1", "t2"])
    A, x = PolynomialRing(R, "x")
    phi1=x^12-A(5*t1^3*t2^4)#v(x) = 1/4+sqrt(2)/3
    phi2=phi1^4-2^9*3^8*t1^11*t2^7*x^18# v(phi1) = 31/8+13/4*sqrt(2)
    F=phi2^3-A(144*3^25*5^14*t1^38*t2^24)# v(phi2) =  38/3+8*sqrt(2) ; 144 = 2^(-858) mod 149

    Rt, z = PolynomialRing(QQ, "z")
    Ka, a = NumberField(z^2-2, "a")

    Phi = [x,phi1,phi2]
    vals = [1//4+a//3,31//8+13//4*a] # v(phi2) not considered here4433

    elt = PhiExp(F,Phi)
    return PhiNewtonPolygon(elt,vals)==[[0,38+24*a],[3,0]]
end

function TestAllCoeffGivenVKt1t2()
    R, (t1, t2) = LaurentPolynomialRing(GF(149), ["t1", "t2"])
    A, x = PolynomialRing(R, "x")
    phi1=x^12-A(5*t1^3*t2^4)#v(x) = 1/4+sqrt(2)/3
    phi2=phi1^4-2^9*3^8*t1^11*t2^7*x^18# v(phi1) = 31/8+13/4*sqrt(2)
    F=phi2^3-A(144*3^25*5^14*t1^38*t2^24)# v(phi2) =  38/3+8*sqrt(2) ; 144 = 2^(-858) mod 149

    Rt, z = PolynomialRing(QQ, "z")
    Ka, a = NumberField(z^2-2, "a")

    Phi = [x,phi1,phi2]
    vals = [1//4+a//3,31//8+13//4*a]

    elt = PhiExp(F,Phi)
    return AllCoeffGivenV(elt[1],vals,38+24*a)==[[GF(149)(114), [38,24], 0, 0]]
end

function TestPhiResidualPolKt1t2()
    p=1523
    Ff0base=FiniteField(p,1)[1]
    # ICI : si j'utilise GF(p), quand je vais bosser avec des extensions de corps au niveau des polynomes residuels,
    # les conversions ne marcheront pas...

    R, (t1, t2) = LaurentPolynomialRing(Ff0base, ["t1", "t2"])
    A, x = PolynomialRing(R, "x")

    Rt, z = PolynomialRing(QQ, "z")
    Ka, a = NumberField(z^2-2, "a")

    phi0=x^2+x+1 # introduit une extension de degre 2 (j)
    phi1=phi0^72+A(1406*t1^6*t2^4) # pente 1/12+1/18*sqrt(2), ramification 36 ; je pense aussi introduire l'extension sqrt(xxx)
    F=phi1^8-113*t1^57*t2^30*phi0^36 # pente 3/2 pour RNP, i.e. 15/2+4*sqrt(2) à la Nart, ramification 4, extention sqrt(yyy) ?

    Phi=[x,phi0,phi1]
    vals=[0,1//12+a//18,15//2+4*a]

    Lambda0=[[Ff0base(1),Ff0base(1)]]
    F0=PhiExp(F,Phi[1:1])
    N0=PhiNewtonPolygon(F0,[])
    e0=1
    R0=PhiResidualPol(F0, N0, vals[1:0], Lambda0, e0)
    Ff0,y = PolynomialRing(Ff0base, "y")
    if R0 != (y^2+y+1)^576 return false end

    F1=PhiExp(F,Phi[1:2])
    N1=PhiNewtonPolygon(F1,vals[1:1])
    e1=36
    Ff0,y = PolynomialRing(GF(p), "y")
    Ff1base,z0 = FiniteField(y^2+y+1,"z0")
    Ff1,y = PolynomialRing(Ff1base, "y")
    Lambda1 = [[z0,z0],Ff1base(1)]

    R1=PhiResidualPol(F1, N1, vals[1:1], Lambda1, e1)
    if R1!=(y^2 + 1406*z0)^8 return false end
    return true
end