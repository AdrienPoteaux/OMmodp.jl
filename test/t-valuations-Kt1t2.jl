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

function Main.OMFacto.:BivCoeff(elt::AbstractAlgebra.Generic.Poly{AbstractAlgebra.Generic.LaurentMPolyWrap{fpFieldElem, fpMPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrapRing{fpFieldElem, fpMPolyRing}}}, i::Int64, n::nf_elem)
    tmp=coeff(elt,i)
    c=representation_matrix(n)[1,:]
    shift=[0,0]
    # There is no coeff function for Multivariate Laurent Polynomials.
    # We need to shift it via .mpoly but I do not know a way to get the shift done by this function...
    # Therefore, I precompute this shift in a not nice way
    for i in exponent_vectors(tmp)
        if i[1] < shift[1] shift[1]=i[1] end
        if i[2] < shift[2] shift[2]=i[2] end
    end
    # Have to make a weird cast here
    return coeff(tmp.mpoly,[Int64(numerator(c[1]))-shift[1],Int64(numerator(c[2]))-shift[2]])
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

## Now we can write an example. 
# This is a big irreducible example I used in our exchange about multivariate Puiseux expansions.
# I assume the reduction modulo 149 is fine.

function TestPhiValKt1t2()
    R, (t1, t2) = LaurentPolynomialRing(GF(149), ["t1", "t2"])
    A, x = PolynomialRing(R, "x")
    phi1=x^12-A(5*t1^3*t2^4)#v(x) = 1/4+sqrt(2)/3
    phi2=phi1^4-2^9*3^8*t1^11*t2^7*x^18# v(phi1) = 31/8+13/4*sqrt(2)
    F=phi2^3-A(144*3^25*5^14*t1^38*t2^24)# v(phi2) =  38/3+8*sqrt(2) ; 144 = 2^(-858) mod 149

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
    return AllCoeffGivenV(elt[1],vals,38+24*a)==[[GF(149)(114), 24*a + 38, 0, 0, 0]]
end
