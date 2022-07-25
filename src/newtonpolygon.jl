function GaussVal(F)
#  In: F in K[[x]][y]
# Out: the Gauss valuation of F (-1 if input is 0)
    tmp = minimum(map(valuation,F.coeffs))
    if tmp == F.parent.base_ring.prec_max
        return -1
    end
    return tmp
end

function PhiVal(elt, vals)
#  In: elt some Phi expansion of F in K[[x]][y], vals the valuation of the key polynomials
# Out: the associated valuation (assumed positive, as we return -1 for the 0 polynomial)
    k=length(vals)
    last=length(elt)
    mini=-1
    while mini == -1
        mini = k>1 ? PhiVal(elt[last],vals[1:k-1]) : GaussVal(elt[last])
        last-=1
        if last==0
            break
        end
    end
    mini += vals[k]*last # This adds 0 if last==0
    for i in 1:last
        tmp = k>1 ? PhiVal(elt[i],vals[1:k-1]) : GaussVal(elt[i])
        if tmp != -1
            tmp+=vals[k]*(i-1)
        end
        if mini>tmp
            mini=tmp
        end
    end
    return mini
end

function NewtonPolygon(exp,vals)
#  In: exp a Phi expansion of some polynomials, vals the valuation of the key polynomials
# Out: the associated generalised Newton polygon
    k=length(vals)
    points=[PhiVal(elt[i],vals[1:k-1]) for i in 1:length(elt[i])]
    # calculer l'enveloppe convexe... voir si ca existe dans Julia ?
end