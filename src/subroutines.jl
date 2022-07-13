function TaylorExp(F, phi)
    #  In: F, phi in K[x]/(x^n)[y]
    # Out: [a_0,...,a_s] s.t. F=a_0+a_1*phi+...+a_s*phi^s
    return _TaylorExp(F,phi)
end

function _TaylorExp(F, phi)
    d=F.length-1
    m=phi.length-1
    k=div(d,m)
    if k<=0
       return F
    end
    k=div(k+1,2) # ceil(k/2)
    tmp=phi^k
    q,r=divrem(F,tmp)
    return [_TaylorExp(r,phi);_TaylorExp(q,phi)]
end