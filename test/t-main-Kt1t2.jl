function TestFirstApproximantsKt1t2()
    p=1523
    F0=GF(p)
    R, (t1, t2) = LaurentPolynomialRing(F0, ["t1", "t2"])
    A, x = PolynomialRing(R, "x")
    phi0=x^2+x+1 # introduit une extension de degre 2 (j)
    phi1=phi0^72+A(1406*t1^6*t2^4) # pente 1/12+1/18*sqrt(2), ramification 36 ; pas d'extension de corps ici -> reductible.
    P=phi1^8-113*t1^57*t2^30*phi0^36 # pente 3/2 pour RNP, i.e. 15/2+4*sqrt(2) à la Nart, ramification 4, extention sqrt(yyy) ?

    res=FirstApproximants(P)
    if res[1]!=false return false end
    res=res[2]
    if length(res)!=2 return false end
    if (degree(res[1][1]) != 1) || (degree(res[2][1]) != 2) return false end
    z = PolynomialRing(QQ, "z")[2]
    Ka, a = NumberField(z^2-2, "a")
    if (res[1][2]!=0) || (res[2][2]!=1//12+a//18) return false end
    y=PolynomialRing(GF(p),"y")[2]
    if res[1][3]!=y^2+1143 return false end
    return true
end

## Work in progress. PhiTruncate is actually very slow, so this would take ages to run.
function TestPhiHenselKt1t2()
    p=1523
    F0=GF(p)
    R, (t1, t2) = LaurentPolynomialRing(F0, ["t1", "t2"])
    A, x = PolynomialRing(R, "x")
    phi0=x^2+x+1 # introduit une extension de degre 2 (j)
    phi1=phi0^72+A(1406*t1^6*t2^4) # pente 1/12+1/18*sqrt(2), ramification 36 ; pas d'extension de corps ici -> reductible.
    P=phi1^8-113*t1^57*t2^30*phi0^36 # pente 3/2 pour RNP, i.e. 15/2+4*sqrt(2) à la Nart, ramification 4, extention sqrt(yyy) ?

    res=FirstApproximants(P) # This should outputs e too. ICI TODO
    if ((res[1]!=false) || (res[3][4].fac.count!=2)) error("Pb FirstApproximants") end
    Phi=res[3][1]
    vals=res[3][2]
    Lambda=res[3][3]
    fac=res[3][4]
    g,h=[i.first^i.second for i in fac]
    v=4*72*vals[2]
    z = PolynomialRing(QQ, "z")[2]
    Ka, a = NumberField(z^2-2, "a")
    return PhiHensel(P,g,h,v,v,Phi,vals,Lambda,36,Ka(12))
end