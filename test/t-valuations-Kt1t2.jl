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
    F0=GF(p)

    R, (t1, t2) = LaurentPolynomialRing(F0, ["t1", "t2"])
    A, x = PolynomialRing(R, "x")

    Rt, z = PolynomialRing(QQ, "z")
    Ka, a = NumberField(z^2-2, "a")

    phi0=x^2+x+1 # introduit une extension de degre 2 (j)
    phi1=phi0^72+A(1406*t1^6*t2^4) # pente 1/12+1/18*sqrt(2), ramification 36 ; je pense aussi introduire l'extension sqrt(xxx)
    P=phi1^8-113*t1^57*t2^30*phi0^36 # pente 3/2 pour RNP, i.e. 15/2+4*sqrt(2) à la Nart, ramification 4, extention sqrt(yyy) ?

    Phi=[x,phi0,phi1]
    vals=[0,1//12+a//18,15//2+4*a]

    Lambda0=[[F0(1),F0(1)]]
    P0=PhiExp(P,Phi[1:1])
    N0=PhiNewtonPolygon(P0,[])
    e0=1
    R0=PhiResidualPol(P0, N0, vals[1:0], Lambda0, e0)
    F0y,y = PolynomialRing(F0, "y")
    if R0 != (y^2+y+1)^576 return false end

    P1=PhiExp(P,Phi[1:2])
    N1=PhiNewtonPolygon(P1,vals[1:1])
    e1=36
    F1,z0 = FiniteField(y^2+y+1,"z0")
    F1y,y = PolynomialRing(F1, "y")
    Lambda1 = [[z0,z0],F1(1)]

    R1=PhiResidualPol(P1, N1, vals[1:1], Lambda1, e1)
    if R1!=(y^2 + 1406*z0)^8 return false end
    return true
end