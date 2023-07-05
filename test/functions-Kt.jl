function Main.OMFacto.:CoeffAndExp(elt::fpRelPowerSeriesRingElem, n::Int64) # ICI GENERIQUE j'aimerais preciser le RingElem en K[[t]]; Comment faire ?
    return [coeff(elt,n),[n]]
end

function Main.OMFacto.:ValueGroup(A::fpRelPowerSeriesRing)
    return [Rational(1)]
end

function Main.OMFacto.:ResidueField(A::fpRelPowerSeriesRing)
    return base_ring(A) # enough ?
end

function Main.OMFacto.:GammaCofactors(Gamma::Vector{Rational{T}},gamma::Rational{T}) where {T}
    return [gamma//Gamma[1]]
end

# In  : F some base field, P an irreducible polynomial above it, s a string to express the new variable.
# Out : The field extension F(P), and an element of this new field that is a root of P.
function Main.OMFacto.IncResField(F::Nemo.fpField,P::fpPolyRingElem,s::String)
    # Input is Fp, so this is easy
    return FiniteField(P,s)
end

# same for Fq
function Main.OMFacto.IncResField(F::fqPolyRepField,P::fqPolyRepPolyRingElem,s::String)
    # if F has degree 1, I am not sure this works.
    p=Int(characteristic(F))
    z=gen(F)
    FF,zz=FiniteField(p,degree(F)*degree(P),s)
    # Should we check that the minimal polynomial of F factorise on FF ?
    y=PolynomialRing(FF,"y")[2]
    tmp=roots(sum([FF(coeff(P,i))*y^i for i in 0:degree(P)])) # not testing if there is a root here...
    return FF, tmp[1]
end
