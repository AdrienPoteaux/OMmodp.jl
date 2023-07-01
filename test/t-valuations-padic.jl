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

function TestPhiResidualPolPadic()
    A = PadicField(3,60) # ICI : y'a un souci avec la precision... on en perd vraiment beaucoup en route
    # (avec precision 50, je trouve pas le bon N0 car le coefficient constant est tronque a precision 15)
    x = PolynomialRing(A,"x")[2]

    P = (x^2-27)^5+2*3^15*x
    Phi=[x,x^2-27]

    vals = [3//2]

    P0 = PhiExp(P,Phi[1:1])
    N0 = PhiNewtonPolygon(P0,[])
    e0=2
    F0 = GF(3)
    y = PolynomialRing(F0,"y")[2]

    Lambda0=[[F0(1)]]
    R0=PhiResidualPol(P0,N0,[],Lambda0,e0)
    if R0!=(y-1)^5 return false end

    z0=F0(1)
    P1 = PhiExp(P,Phi)
    N1 = PhiNewtonPolygon(P1,vals)
    e1=5
    F1 = F0
    Lambda1=[[F1(1)],F1(1)]

    R1=PhiResidualPol(P1,N1,vals,Lambda1,e1)

    if R1!=(y-1) return false end

    return true
end