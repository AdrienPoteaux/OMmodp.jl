## Je prends l'exemple des slides pour rester simple

function Main.OMFacto.:CoeffAndExp(elt::padic, n::Int64)
    aq=lift(QQ,elt)
    p=elt.parent.p
    return [floor(mod(aq,p^(n+1))//p^n),[n]]
end

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
    return AllCoeffGivenV(elt[1],vals,33//2) == [[2,[15],1]]
end