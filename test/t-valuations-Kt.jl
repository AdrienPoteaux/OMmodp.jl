Main.OMFacto.eval(:(
    function ConstInv(A::Generic.RelSeries, n::Int64)
        return base_ring(A)(1)//n
    end
    )
)

Main.OMFacto.eval(:(
    function BivCoeff(elt::Generic.Poly{fpRelPowerSeriesRingElem}, i::Int64, n::Int64) # ICI GENERIQUE j'aimerais preciser le RingElem en K[[t]]; Comment faire ?
        return coeff(coeff(elt,i),n)
    end
    )
)

function TestPhiValKt()
    A, t = PowerSeriesRing(GF(211), 33, "t")
    L, x = PolynomialRing(A,"x")

    F=x^8+(t^10+203*t^2)*x^7+(28*t^4+201*t^3)*x^6+(143*t^6+60*t^5)*x^5+(209*t^9+130*t^8+61*t^7+33*t^6)*x^4+(8*t^11+35*t^10+37*t^9+79*t^8)*x^3+(190*t^12+128*t^11+198*t^10+171*t^9)*x^2+(7*t^12+80*t^11)*x+16*t^12
    Phi=[x+132*t^10+210*t^2,x^2+(53*t^10+209*t^2)*x+107*t^12+t^4+210*t^3,x^4+(106*t^10+207*t^2)*x^3+(2*t^12+6*t^4+206*t^3)*x^2+(5*t^14+108*t^13+201*t^6+10*t^5)*x+15*t^15+210*t^9+7*t^8+206*t^7+4*t^6]
    vals=[3//2,15//4,15//2] # last one not augmented ?
    elt=PhiExp(F,Phi)
    return PhiVal(elt,vals) == 13
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
    return AllCoeffGivenV(elt[1],vals,13) == [[GF(211)(40),13,0,0,0]]
end