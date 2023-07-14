Main.OMFacto.eval(:(
    function ConstInv(A::GFPRelSeriesRing, n::Int64)
        return base_ring(A)(1)//n
    end
    )
)

function TestAppRootKt()
    A, x = PowerSeriesRing(GF(211), 33, "x")
    L, y = PolynomialRing(A,"y")
    F=y^4+x^3*y^3+209*x^3*y^2+x^6

    phi=AppRoot(F,1)
    if (F-phi != 0)
        return false
    end

    phi=AppRoot(F,4)

    if (degree(F-phi^4)>2)
        return false
    end

    phi=AppRoot(F,2)

    if (degree(F-phi^2)>2)
        return false
    end

    return true
end

function TestTaylorExpKt()
    A, x = PowerSeriesRing(GF(211), 33, "x")
    L,y=PolynomialRing(A,"y")
    phi = y^2+(100*x^3+x^2+4)*y+2*x^2+105*x+98
    F = (53*x^12+102*x^9+65*x^6)*y^4+(160*x^9+146*x^6)*y^3+(70*x^6+204*x^3)*y^2+(109*x^3+32*x^2+5*x)*y+3*x+11
    tmp=TaylorExp(F,phi)
    somme=0
    for i in eachindex(tmp)
        somme=somme+phi^(i-1)*tmp[i]
    end
    
    if somme != F
        return false
    end

    phi=y^2+(53*x^10+209*x^2)*y+107*x^12+x^4+103*x^3
    F=y^8+(x^10+203*x^2)*y^7+(28*x^4+207*x^3)*y^6+(135*x^6+24*x^5)*y^5+(194*x^9+170*x^8+151*x^7+6*x^6)*y^4+(68*x^11+166*x^10+120*x^9+187*x^8)*y^3+(109*x^13+115*x^12+31*x^11+36*x^10+207*x^9)*y^2+(148*x^15+118*x^14+144*x^13+167*x^12+8*x^11)*y+150*x^15+26*x^14+207*x^13+x^12

    tmp=TaylorExp(F,phi)
    somme=0
    for i in eachindex(tmp)
        somme=somme+phi^(i-1)*tmp[i]
    end
    
    if somme != F
        return false
    end

    return true
end

function TestPhiExpKt()
    A, x = PowerSeriesRing(GF(211), 33, "x")
    L, y = PolynomialRing(A,"y")
    F=y^8+(x^10+203*x^2)*y^7+(28*x^4+201*x^3)*y^6+(143*x^6+60*x^5)*y^5+(209*x^9+130*x^8+61*x^7+33*x^6)*y^4+(8*x^11+35*x^10+37*x^9+79*x^8)*y^3+(190*x^12+128*x^11+198*x^10+171*x^9)*y^2+(7*x^12+80*x^11)*y+16*x^12
    Phi=[y+132*x^10+210*x^2,y^2+(53*x^10+209*x^2)*y+107*x^12+x^4+103*x^3,y^4+(106*x^10+207*x^2)*y^3+(2*x^12+6*x^4+206*x^3)*y^2+(5*x^14+108*x^13+201*x^6+10*x^5)*y+15*x^15+210*x^9+7*x^8+206*x^7+4*x^6]
    res=PhiExp(F,Phi)
    if PhiEval(res,Phi) != F
        return false
    end

    F=y^8+(x^10+203*x^2)*y^7+(28*x^4+207*x^3)*y^6+(135*x^6+24*x^5)*y^5+(194*x^9+170*x^8+151*x^7+6*x^6)*y^4+(68*x^11+166*x^10+120*x^9+187*x^8)*y^3+(109*x^13+115*x^12+31*x^11+36*x^10+207*x^9)*y^2+(148*x^15+118*x^14+144*x^13+167*x^12+8*x^11)*y+150*x^15+26*x^14+207*x^13+x^12
    Phi=[y+132*x^10+210*x^2,y^2+(53*x^10+209*x^2)*y+107*x^12+x^4+210*x^3]
    res=PhiExp(F,Phi)
    if PhiEval(res,Phi) != F
        return false
    end

    return true
end