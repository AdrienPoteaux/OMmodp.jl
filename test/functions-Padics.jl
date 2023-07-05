function Main.OMFacto.:CoeffAndExp(elt::padic, n::Int64)
    aq=lift(QQ,elt)
    p=elt.parent.p
    return [floor(mod(aq,p^(n+1))//p^n),[n]]
end