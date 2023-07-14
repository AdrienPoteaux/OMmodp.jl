function BaseGenerators(A)
    error("BaseGenerators must be defined for our base ring. It is supposed to provide generators of the base ring")
end

function ResidueField(A)
    error("ResidueField must be defined ; input type: ",typeof(A))
end

# Je n'arrive pas a voir comment faire correctement les choses de maniere generique pour ce point
# Neanmoins, ca doit etre possible.
function GammaCofactors(Gamma::Vector,gamma)
    error("GammaCofactors must be defined. Given an element gamma and a basis of Gamma, this outputs a vector v such that gamma=Gamma.v (scalar product)")
end

function IncResField(F,R,s::String)
    error("IncResField must be defined. Input types :",typeof(F)," ; ",typeof(R)," ; ",typeof(s))
end

# An xgcd that ensures 0 <= lb < b/gcd(a,b)
# (probably not mandatory, but ensures similar results at each call)
function mygcdx(a,b)
    g,l,lb=gcdx(a,b)
    while l<0
        l+=b
        lb-=a
    end
    while l>=(b//g)
        l-=b
        lb+=a
    end
    return g,l,lb
end

# Computes the slope and associated data to the edge Delta.
#  In : Delta the slope (a couple of points), OldGamma the base of the current group value
# Out : A tuple containing
#       * gamma the new augmented valuation (minus the slope)
#       * Gamma a base of the new group value
#       * e the ramification index associated to the slope
#       * Lb the product of the l'_j (cf NaOl21)
#       * L the vector of the L_j (cf NaOl21)
function OneSlope(Delta::Vector, OldGamma::Vector)
    gamma=-(Delta[2][2]-Delta[1][2])//(Delta[2][1]-Delta[1][1])
    v=GammaCofactors(OldGamma,gamma)
    k=length(OldGamma)
# on va declarer plein de vecteurs pour se simplifier la vie.
# Possible que ce soit optimisable...
    h=[numerator(v[j]) for j in eachindex(v)]
    e=[denominator(v[j]) for j in eachindex(v)]
    eb=[1 for j in 1:k]
    d=[1 for j in 1:k]
    l=[0 for j in 1:k]
    lb=[0 for j in 1:k]
    prode=1
    prodd=1
    d[1],l[1],lb[1]=mygcdx(h[1]*eb[1],e[1])
    for j in 2:k
        prodd=prodd*d[j-1]
        prode=prode*e[j-1]
        eb[j]=numerator(prode//prodd)
        d[j],l[j],lb[j]=mygcdx(h[j]*eb[j],e[j])
    end
    prodd=prodd*d[k]
    prode=prode*e[k]
    
    # have to deal with j=k separately to avoid sum([]) effect.
    Gamma=[[d[j]//e[j]*OldGamma[j]+sum([l[j]*eb[j]*h[m]//e[m]*OldGamma[m] for m in j+1:k]) for j in 1:k-1];d[k]//e[k]*OldGamma[k]]
    L=[[l[j]*prod([lb[m] for m in j+1:k]) for j in 1:k-1];l[k]]
    Lb=prod(lb)
    return gamma,Gamma,numerator(prode//prodd),Lb,L
end

"""
    FirstApproximants(P::Generic.Poly{T}) where {T}

    Test if the polynomial P is irreducible over its base ring.
    If it is, the function outputs true, together with the list of successive
    (key polynomial, associated polynomial, irreducible factor of the residual polynomial)

    If it is not, the same list is provided (up to the point where we prove non irreducibility)
    A third element is returned : [Phi,vals,Lambda,NP or RP] where
     * Phi is the list of key Polynomials
     * vals is the list associated augmented valuations
     * Lambda is the current list of correcting terms (for residual polynomial computation)
     * Last parameter is either (the size of the last tuple in the 2nd output parameter tells us which one)
        ** the last computed Newton polygon (with several slopes)
        ** the factorisation of the last computed residual polynomial (with at least two coprime factors)
    This list will be provided to the Hensel algorithm to start the factorisation of P

    IMPORTANT : several functions are supposed to exist for our base ring
     * A valuation function (elements of the value groupe must be provided with a way to compare them ; one must also have a numerator function for them).
     * A function CoeffAndExp.
     * BaseGenerators a function that provides the generators of the base ring.
     * ResidueField that provides the residue field associated to the base ring.
     * GammaCofactors that given a base for some ZZ-group Gamma and an element gamma outputs a vector v (of rationals) such that gamma=Gamma.v (scalar product)
     * IncResField a function that given a field Ff and an irreducible polynomial Pf over it, compute the associated field extension and a root of Pf in this field
"""
# ICI : we need to add an optional parameter that provide already computed things for recursive calls after some factorisation step.
# Le fait de donner approximants en sortie quand c'est pas irreductible n'est utile que pour le reutiliser en cas de futur appel recursif ?
function FirstApproximants(P::Generic.Poly{T}) where {T}
    N=degree(P)
    A=base_ring(P)
    FirstGamma=map(valuation,BaseGenerators(A)) # we need to keep that one for later use.
    k=length(FirstGamma)
    Gamma=FirstGamma
    vals=[Gamma[1]][1:0] # so that type is always correct
    Phi=[P][1:0] # same
    F=ResidueField(A)
    Lambda=[[F(1) for i in 1:k]] # This is the list of R_0(\tau_{-1,j})
#    z=F(1) # inutile ?
    approximants=[] # will be the list of computed approximants ; list of phi,gamma,R
    while N>1
        Phi=[Phi;AppRoot(P,N)]
        dvt=PhiExp(P,Phi)
        NP=PhiNewtonPolygon(dvt,vals)
        approximants=[approximants;[last(Phi)]]
        if length(NP)>2 return [false,approximants,[Phi,vals,Lambda,NP]] end
        gamma,Gamma,e,Lb,L=OneSlope(NP, Gamma)
        approximants[length(approximants)]=[approximants[length(approximants)];gamma]
        RP=PhiResidualPol(dvt, NP, vals, Lambda, e)
        vals=[vals;gamma]
        fac=factor(RP)
        if fac.fac.count > 1 return [false,approximants,[Phi,vals,Lambda,fac]] end
        tmp=[i for i in fac.fac][1]
        R=tmp.first
        N=tmp.second
        approximants[length(approximants)]=[approximants[length(approximants)];R]
        if degree(R) > 1
            F,z=IncResField(F,R,string("z",length(vals)-1))
            # converting correcting elements in the new residue field.
            Lambda=[[[F(i) for i in Lambda[1]]];[F(Lambda[i]) for i in 2:length(Lambda)]]
        else 
            z=-coeff(R,0) # assuming R is monic here
        end
        # Now multiplying elements of Lambda by the appropriate powers of z.
        for j in 1:k
            tmp=GammaCofactors(Gamma,FirstGamma[j])
            # la sortie dans tmp c'est des rationnels
            # mais ici on sait que ce sont des entiers (car on decompose un element de Gamma_i dans une base de Gamma_i)
            # on passe donc par numerator
            Lambda[1][j]=Lambda[1][j]*z^(-Int(numerator(sum([tmp[i]*L[i] for i in 1:k]))))
        end
        for i in 2:length(Lambda)
            tmp=GammaCofactors(Gamma,vals[i-1])
            Lambda[i]=Lambda[i]*z^(-numerator(sum([tmp[i]*L[i] for i in 1:k])))
        end
        Lambda=[Lambda;z^Lb]
    end
    return [true,approximants]
end