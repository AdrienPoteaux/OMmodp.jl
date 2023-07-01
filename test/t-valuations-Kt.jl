function Main.OMFacto.:CoeffAndExp(elt::fpRelPowerSeriesRingElem, n::Int64) # ICI GENERIQUE j'aimerais preciser le RingElem en K[[t]]; Comment faire ?
    return [coeff(elt,n),[n]]
end

# same for Fq
function Main.OMFacto.:CoeffAndExp(elt::fqPolyRepRelPowerSeriesRingElem, n::Int64)
    return [coeff(elt,n),[n]]
end

function TestPhiValKt()
    A, t = PowerSeriesRing(GF(211), 33, "t")
    L, x = PolynomialRing(A,"x")

    F=x^8+(t^10+203*t^2)*x^7+(28*t^4+207*t^3)*x^6+(135*t^6+24*t^5)*x^5+(194*t^9+170*t^8+151*t^7+6*t^6)*x^4+(68*t^11+166*t^10+120*t^9+187*t^8)*x^3+(109*t^13+115*t^12+31*t^11+36*t^10+207*t^9)*x^2+(148*t^15+118*t^14+144*t^13+167*t^12+8*t^11)*x+150*t^15+26*t^14+207*t^13+t^12
    Phi=[x+132*t^10+210*t^2,x^2+(53*t^10+209*t^2)*x+107*t^12+t^4+210*t^3,x^4+(106*t^10+207*t^2)*x^3+(2*t^12+6*t^4+206*t^3)*x^2+(5*t^14+108*t^13+201*t^6+10*t^5)*x+15*t^15+210*t^9+7*t^8+206*t^7+4*t^6]
    vals=[3//2,15//4,6] # la 3eme correspond à la plus petite pente de N2(F) = (0,51/4),(1,6),(2,0)
    elt=PhiExp(F,Phi)
    return PhiVal(elt,vals) == 12
end

function TestPhiNewtonPolygonKt()
    A, t = PowerSeriesRing(GF(211), 33, "t")
    L ,x = PolynomialRing(A,"x")

    F=x^8+(t^10+203*t^2)*x^7+(28*t^4+207*t^3)*x^6+(135*t^6+24*t^5)*x^5+(194*t^9+170*t^8+151*t^7+6*t^6)*x^4+(68*t^11+166*t^10+120*t^9+187*t^8)*x^3+(109*t^13+115*t^12+31*t^11+36*t^10+207*t^9)*x^2+(148*t^15+118*t^14+144*t^13+167*t^12+8*t^11)*x+150*t^15+26*t^14+207*t^13+t^12
    Phi=[x+132*t^10+210*t^2,x^2+(53*t^10+209*t^2)*x+107*t^12+t^4+210*t^3]
    vals=[3//2]
    elt=PhiExp(F,Phi)
    return PhiNewtonPolygon(elt,vals)==[[0,15],[4,0]]
end

function TestAllCoeffGivenVKt()
    A, t = PowerSeriesRing(GF(211), 33, "t")
    L, x = PolynomialRing(A,"x")

    F=x^8+(t^10+203*t^2)*x^7+(28*t^4+201*t^3)*x^6+(143*t^6+60*t^5)*x^5+(209*t^9+130*t^8+61*t^7+33*t^6)*x^4+(8*t^11+35*t^10+37*t^9+79*t^8)*x^3+(190*t^12+128*t^11+198*t^10+171*t^9)*x^2+(7*t^12+80*t^11)*x+16*t^12
    Phi=[x+132*t^10+210*t^2,x^2+(53*t^10+209*t^2)*x+107*t^12+t^4+103*t^3,x^4+(106*t^10+207*t^2)*x^3+(2*t^12+6*t^4+206*t^3)*x^2+(5*t^14+108*t^13+201*t^6+10*t^5)*x+15*t^15+210*t^9+7*t^8+206*t^7+4*t^6]
    elt=PhiExp(F,Phi)   
    vals=[3//2,15//4]
    return AllCoeffGivenV(elt[1],vals,13) == [[GF(211)(40),[13],0,0]]
end

function TestPhiResidualPolKt()
    p = 211
    F0 = FiniteField(p,1)[1]
    A, t = PowerSeriesRing(F0, 33, "t")
    L, x = PolynomialRing(A,"x")

    phi0=x^2+x+3
    P=phi0^4-3*t^2
    Phi=[x,phi0]

    P0=PhiExp(P,Phi[1:1])
    N0=PhiNewtonPolygon(P0,[])
    e0=1
    vals=[0]
    Lambda0=[[F0(1),F0(1)]]
    R0=PhiResidualPol(P0,N0,[],Lambda0,e0)

    y=PolynomialRing(F0,"y")[2]
    if R0!=(y^2+y+3)^4 return false end

    y=PolynomialRing(GF(p),"y")[2] # ICI faudrait voir comment gerer les extensions de Fq dnas le futur
    F1,z0=FiniteField(y^2+y+3,"z0")
    P1=PhiExp(P,Phi)
    N1=PhiNewtonPolygon(P1,vals)
    e1=2
    # applying formulae, we get :
    Lambda1=[[F1(1)],z0]
    R1=PhiResidualPol(P1,N1,vals,Lambda1,e1)

    y=PolynomialRing(F1,"y")[2] # ICI faudrait voir comment gerer les extensions de Fq dnas le futur
    if R1!=y^2+208 return false end

    return true
end