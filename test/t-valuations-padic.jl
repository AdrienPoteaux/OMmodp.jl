## Je prends l'exemple des slides pour rester simple

Main.OMFacto.eval(:(
    function ConstInv(A::GFPRelSeriesRing, n::Int64)
        return base_ring(A)(1)//n
    end
    )
)


Main.OMFacto.eval(:(
    function BivCoeff(elt::Generic.Poly{padic}, i::Int64, n::Int64)
        a=coeff(elt,i)
        aq=lift(QQ,a)
        p=a.parent.p
        return floor(mod(aq,p^(n+1))//p^n)
    end
    )
)

function TestPhiValPadic()
    A    = PadicField(3, 50)
    L, x = PolynomialRing(A,"x")

    F = (x^2-27)^5+3^15*x## le x^10 part dans le warp ???
    Phi = [x,x^2-27]
    vals=[3//2,33//10]
    elt=PhiExp(F,Phi)
    return PhiVal(elt,vals) == 33//2
end

function TestPhiNewtonPolygonPadic()
    A    = PadicField(3, 50)
    L, x = PolynomialRing(A,"x")

    F = (x^2-27)^5+3^15*x
    Phi = [x,x^2-27]
    vals=[3//2]
    elt=PhiExp(F,Phi)
    return PhiNewtonPolygon(elt,vals)==[[0,33//2],[5,0]]
end

function TestAllCoeffGivenVPadic()
    A    = PadicField(3, 50)
    L, x = PolynomialRing(A,"x")

    F = (x^2-27)^5+2*3^15*x # just modified the last coefficient here
    Phi = [x,x^2-27]
    vals=[3//2]
    elt=PhiExp(F,Phi)
    return AllCoeffGivenV(elt[1],vals,33//2) == [[2,15,0,1]]
end