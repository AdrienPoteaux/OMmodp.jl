using OMmodp
using Nemo

function TestTaylorExp()
    res=true
    R, x = PolynomialRing(GF(211), "x")
    T=ResidueRing(R, x^33) # needed ?
    L,y=PolynomialRing(T,"y")
    phi = y^2+(100*x^3+x^2+4)*y+2*x^2+105*x+98
    F = (53*x^12+102*x^9+65*x^6)*y^4+(160*x^9+146*x^6)*y^3+(70*x^6+204*x^3)*y^2+(109*x^3+32*x^2+5*x)*y+3*x+11
    tmp=TaylorExp(F,phi)
    somme=0
    for i in 1:length(tmp)
        somme=somme+phi^(i-1)*tmp[i]
    end
    
    if somme != F
        res=false
    end

    phi=y^2+(53*x^10+209*x^2)*y+107*x^12+x^4+103*x^3
    F=y^8+(x^10+203*x^2)*y^7+(28*x^4+207*x^3)*y^6+(135*x^6+24*x^5)*y^5+(194*x^9+170*x^8+151*x^7+6*x^6)*y^4+(68*x^11+166*x^10+120*x^9+187*x^8)*y^3+(109*x^13+115*x^12+31*x^11+36*x^10+207*x^9)*y^2+(148*x^15+118*x^14+144*x^13+167*x^12+8*x^11)*y+150*x^15+26*x^14+207*x^13+x^12

    tmp=TaylorExp(F,phi)
    somme=0
    for i in 1:length(tmp)
        somme=somme+phi^(i-1)*tmp[i]
    end
    
    if somme != F
        res=false
    end

    return res
end