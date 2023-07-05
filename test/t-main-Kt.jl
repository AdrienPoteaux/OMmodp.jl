function TestFirstApproximantsKt()
    p=263
    K,t=PowerSeriesRing(GF(p),50,"t")
    x=PolynomialRing(K,"x")[2]
    phi0=x^2+x+1 # f=2 [y^2+y+1]
    phi1=phi0^3-t^2 # mu1(phi0)=2/3 (e=3, f=1 [y-1])
    phi2=phi1^4+2*phi0^2*t^3*phi1^2+(x+132)*t^8*phi0 # mu2(phi1)=13/6 (e=2, f=2)
    P=phi2^5+t^44 # mu3(phi2)=44/5 (e=5, f=1)
    res=FirstApproximants(P)
    if res[1]!=true return false end
    if [degree(i[1]) for i in res[2]]!=[1,2,6,24] return false end
    if [degree(i[3]) for i in res[2]]!=[2,1,2,1] return false end
    if [i[2] for i in res[2]]!=[0,2//3,13//6,44//5] return false end
    return true
end