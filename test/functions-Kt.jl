function Main.OMFacto.:CoeffAndExp(elt::fpRelPowerSeriesRingElem, n::Int64) # ICI GENERIQUE j'aimerais preciser le RingElem en K[[t]]; Comment faire ?
    return [coeff(elt,n),[n]]
end