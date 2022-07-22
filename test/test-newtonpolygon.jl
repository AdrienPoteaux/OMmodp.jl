function TestPhiVal()
    A, x = PowerSeriesRing(GF(211), 33, "x")
    L,y=PolynomialRing(A,"y")

    F=y^8+(x^10+203*x^2)*y^7+(28*x^4+201*x^3)*y^6+(143*x^6+60*x^5)*y^5+(209*x^9+130*x^8+61*x^7+33*x^6)*y^4+(8*x^11+35*x^10+37*x^9+79*x^8)*y^3+(190*x^12+128*x^11+198*x^10+171*x^9)*y^2+(7*x^12+80*x^11)*y+16*x^12
    Phi=[y+132*x^10+210*x^2,y^2+(53*x^10+209*x^2)*y+107*x^12+x^4+210*x^3,y^4+(106*x^10+207*x^2)*y^3+(2*x^12+6*x^4+206*x^3)*y^2+(5*x^14+108*x^13+201*x^6+10*x^5)*y+15*x^15+210*x^9+7*x^8+206*x^7+4*x^6]
    vals=[3//2,15//4,15//2] # last one not augmented ?
    elt=PhiExp(F,Phi)
    return PhiVal(elt,vals) == 13
end