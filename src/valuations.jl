function _turn_left(A,B,C)
#  In: 3 points
# Out: true if one turn left by going from A to C through B
  return ((B[1]-A[1])*(C[2]-A[2])) > ((B[2]-A[2])*(C[1]-A[1]))
end

# not great but...
function _get_2nd_elt(s)
#  In: a stack with at least two elements
# Out: the second element of the stack
  tmp = pop!(s)
  res = first(s)
  push!(s,tmp)
  return res
end

function LowerConvexHull(L) # based on Graham algorithm
#  In: a list L of points (i,j) (we assume that no two points are the same)
# Out: the lower convex hull of this list, assuming L[1] being the origin (with the smallest abscissa)
  # sort -> lt=XXX pour preciser la fonction de tri (qui doit renvoyer true or false ; checker le detail)
  if length(L) <= 2
    return L
  end
  stop = L[length(L)][1] # if the input has increasing abscissae, we store the end to compute the *lower* convex hull only.
  # comparison function
  function cmp(A,B)
    O = L[1]
    turn = (A[1]-O[1])*(B[2]-O[2]) - (A[2]-O[2])*(B[1]-O[1])
    if turn == 0 # A, 0 ,B on the same line.
      return A[1] < B[1]
    end
    return turn > 0 # turn > 0 means OA < OB
  end
  sorted = [[L[1]] ; sort(L[2:length(L)],lt=cmp)]
  s=Stack{parent(L[1])}()
  push!(s,sorted[1])
  push!(s,sorted[2])
  i=3
  while i<=length(L)
    if length(s) == 1
      push!(s,sorted[i])
      i+=1
      continue
    end
    curr = first(s)
    if (curr[1] == stop) break end # we found the lower convex hull
    prec = _get_2nd_elt(s)
    if _turn_left(prec,curr,L[i])
      push!(s,L[i])
      i+=1
    else
      pop!(s)
    end
  end
  # storing the stack in a tab
  res = [pop!(s)]
  while !isempty(s)
    res = [[pop!(s)] ; res]
  end
  return res
end

#= """
    GaussVal(F::Generic.Poly)

  Compute the Gauss valuation of F
  The base ring is supposed to possess a function valuation.
"""
function GaussVal(F::Generic.Poly)
#  In: F in A[y]
# Out: the Gauss valuation of F (-1 if input is 0)
  if F == 0
    return -1
  end
  d=degree(F)
  nzcoeffs=[] # had some examples where F.coeffs was returning a huge list of 0 coefficients, so proceeding that way
  for i in 0:d
    if coeff(F,i) != 0
      nzcoeffs=[nzcoeffs;coeff(F,i)]
    end
  end
  return minimum(map(valuation,nzcoeffs))
end
 =#

 """
    PhiVal(elt, vals::Vector)

Compute the valuation of elt

# arguments
- elt either a table representing some Phi-adic expansion (i.e. a table of table of ... of elements of A), as given by PhiExp
      either an element of A (but then vals must be the empty list)
- vals the (augmented) valuations of the Phi[i]
"""
function PhiVal(elt, vals::Vector)
#  In: elt some Phi expansion of F in A[x], vals the valuation of the key polynomials
# Out: the associated valuation (assumed positive, as we return -1 for the 0 polynomial)
  k=length(vals)
  if k == 0
    # elt is an element of A
    return valuation(elt)
  end
  last=length(elt)
  mini=-1
  while mini == -1
    mini = PhiVal(elt[last],vals[1:k-1])
    last-=1
    if last==0
      break
    end
  end
  mini += vals[k]*last # This adds 0 if last==0
  for i in 1:last
    tmp = PhiVal(elt[i],vals[1:k-1])
    if tmp != -1
      tmp+=vals[k]*(i-1)
      if mini>tmp
        mini=tmp
      end
    end
  end
  return mini
end

"""
    PhiNewtonPolygon(elt::Vector, vals::Vector)

Generalised Newton polygon computation ("Nart" version: valuation of phi_k is not used)

# arguments
- elt a table representing some Phi-adic expansion (i.e. a table of table of ... of polynomials), as given by PhiExp
- vals the (augmented) valuations of the Phi[i]
"""
function PhiNewtonPolygon(elt::Vector, vals::Vector)
#  In: elt a Phi expansion of some polynomials, vals the valuation of the key polynomials
# Out: the associated generalised Newton polygon
  # Remark : we might want to change the implementation to get the "classical" Newton polygon
  valuations=[PhiVal(elt[i],vals) for i in eachindex(elt)] # Nart version for valuations
  first=1 # will be the first point of the list we need to consider
  while (valuations[first] == -1) && (first <= length(valuations))
    first+=1
  end
  if first == length(valuations)+1
    return [] # elt corresponds to some zero polynomial
  end
  l = first # will be leftmost point with smallest ordinate
  v = valuations[first];# will be the valuation of the polynomial
  for i in first+1:length(valuations)
    if (valuations[i] != -1) && (valuations[i] < v)
      l = i
      v = valuations[i]
    end
  end
  last = length(valuations) # will be the last point of the list we need to consider
  while (valuations[last] == -1) # cannot go infinite (we would have already returned [] as a Newton polygon)
    last-=1
  end
  r = last # will be the rightmost point with smallest ordinate
  while (valuations[r] != v)
    r-=1
  end
  K = parent(valuations[1])
  # Now storing the points we have to consider for the lowest convex hull
  left_points = [[K(first-1),valuations[first]]] # will contain the negative slopes
  if first != l
    tmp = (valuations[l]-valuations[first])//(l-first)
    a = -numerator(tmp)
    b = denominator(tmp)
    c = a*(first-1) + b*valuations[first]
    # the line from [first-1,valuations[first]] to [l-1,valuations[l]] is a*i+b*j=c
    for i in first+1:l-1
      if (valuations[i] != -1) && (a*(i-1)+b*valuations[i]<c) # this point can be in the lowest convex hull (-1 means a 0 coefficient, so is not considered).
        left_points = [left_points;[[K(i-1),valuations[i]]]]
      end
    end
    left_points = [left_points;[[K(l-1),valuations[l]]]]
  end
  points = LowerConvexHull(left_points)
  if l!=r
    points = [points;[[K(r-1),valuations[r]]]]
  end
  if r != last
    right_points=[[K(r-1),valuations[r]]]
    tmp = (valuations[r]-valuations[last])//(r-last)
    a = -numerator(tmp)
    b = denominator(tmp)
    c = a*(r-1) + b*valuations[r]
    # the line from [r-1,valuations[r]] to [last-1,valuations[last]] is a*i+b*j=c
    for i in r+1:last-1
      if (valuations[i] != -1) && (a*(i-1)+b*valuations[i]<c) # this point can be in the lowest convex hull.
        right_points=[right_points;[[K(i-1),valuations[i]]]]
      end
    end
    right_points=[right_points;[[K(last-1),valuations[last]]]]
    tmp = LowerConvexHull(right_points)
    points=[points;tmp[2:length(tmp)]]
  end
  return points
end

function CoeffAndExp(elt,n)
  error("CoeffAndExp must be defined for our base ring ; elt type is: ", typeof(elt),n,typeof(n))
end

# More an intern function than anything else
# ICI : ajouter un booleen qui renvoie tous les exposants possibles, peu importe le coefficient ?
function AllCoeffGivenV(elt, vals::Vector, v)
#  In: elt a Phi expansion of some polynomial (either a list of list of ... elements of A, either an element of A)
#      vals the valuation of the key polynomials and v some targeted valuation
# Out: all "coefficients" having valuation v, taking into account all phi.
# One "coefficient" is a tuple [c,[d_{-1,1},...,d_{-1,k}],d_0,...,d_i] such that
# c*pi_{-1,1}^{d_{-1,1}}*...*pi_{-1,k}^{d_{-1,k}}*phi_1^{d_1}*...*phi_i^{d_i} is the "element" of F with targeted valuation.
  k = length(vals)
  res=[]
  if k>0
    for i in eachindex(elt)
      tmp = AllCoeffGivenV(elt[i], vals[1:k-1], v-(i-1)*vals[k])
      for j in eachindex(tmp)
        tmp[j] = [tmp[j] ; i-1]
      end
      if length(tmp)>0 # necessaire ?
        res = [res ; tmp]
      end
    end
    return res
  end
  # k == 0 ; we have an element of A
  if denominator(v)!=1 # no element of A has rational valuation
    return []
  end
  tmp=CoeffAndExp(elt,numerator(v)) # outputs the coefficients and a decomposition of numerator(v) in the fg group \Gamma_{-1} (as a list).
  if tmp[1] != 0
    res=[res;[tmp]]
  end
  return res
end

"""
    PhiResidualPol(elt::Vector, Delta::Vector, vals::Vector, Lambda::Vector, e::Int)

Compute the residual polynomial of elt according to the edge Delta.

vals provide the list of augmented valuations
Lambda the list of "correcting terms" : [[R(pi_{-1,1}), ... , R_(pi_{-1,k})],R(phi_0),...,R(phi_{i-1})]
e the ramification index of the slope.
ResField : the base field of the residual polynomial

IMPORTANT : we assume that a function CoeffAndExp (T, n) exists (the user must define it)
(it is used in the internal function AllCoeffGivenV)
"""
# ICI : je passe par une pente ici, et pas uniquement par une valuation... que veut-on faire in fine ?
# (valuation definie par la dite pente bien sur)
## ICI : we should also code the special case when Delta has a single point.
function PhiResidualPol(elt::Vector, Delta::Vector, vals::Vector, Lambda::Vector, e)
  slope = -(Delta[2][2]-Delta[1][2])//(Delta[2][1]-Delta[1][1])
  ResField=parent(Lambda[1][1])
  Ff,y = PolynomialRing(ResField, "y")
  res=Ff(0)
  for i in numerator(Delta[1][1]):e:numerator(Delta[2][1])
    j = numerator((i-Delta[1][1])//e) # e divides i-Delta[1][1]
    tmp = AllCoeffGivenV(elt[(Int)(i+1)],vals,Delta[1][2]-(j*e*slope))
    for c in tmp
      if length(c) == 2
        res = res + ResField(c[1])*prod([Lambda[1][l]^c[2][l] for l in eachindex(c[2])])*y^(Int(numerator((i-Delta[1][1])//e)))
      else
        res = res + ResField(c[1])*prod([Lambda[1][l]^c[2][l] for l in eachindex(c[2])])*prod([Lambda[l-1]^c[l] for l in 3:length(c)])*y^(Int(numerator((i-Delta[1][1])//e)))
# JULIA : getting a warning that I should use eachindex for the loop 3:length(c)... but I dunno how to avoid indices 1 and 2 with eachindex
      end
    end
  end
  return res
end

# Brut force version, not optimised at all.
function AllMonomialsGivenV(maxdegs::Vector,vals::Vector,Gamma::Vector,v)
  #  In: maxdegs the maximum exponent for each key polynomial, vals the valuation of the key polynomials
  #      Gamma the initial value group and v some targeted valuation
  # Out: all tuple of exponents giving the valuation v, i.e. [d_{-1,1},...,d_{-1,k}],d_0,...,d_i] s.t.
  # val(pi_{-1,1}^{d_{-1,1}}*...*pi_{-1,k}^{d_{-1,k}}*phi_0^{d_0}*...*phi_i^{d_i} = v and d_i<maxdegs[i+1] for i>=0
  k = length(vals)
  if k>0
    res=[]
    for i in 0:maxdegs[k]
      tmp = AllMonomialsGivenV(maxdegs[1:k],vals[1:k-1], Gamma, v-i*vals[k])
      # We add i to each result. Have to make all at once or type conversion are done.
      tmp = [[tmp[j] ; i] for j in eachindex(tmp)]
      res = [res ; tmp]
    end
    return res
  end
  # k == 0 ; we have an element of A
  if denominator(v)!=1 return [] end # no element of A has rational valuation
  return [[map(numerator,GammaCofactors(Gamma,v))]] # GammaCofactors outputs rational numbers
end

"""
    Representant(R,v,Phi::Vector,vals::Vector,Lambda::Vector) where {T}

In : R a residual polynomial (we assume the modified one is given - see PoWe22, Definition 4), v some valuation
     Phi the list of key polynomials, vals the list of their valuation, and Lambda the list of their residual polynomials
Output : A polynomial that has the prescribed residual polynomial and valuation.
"""
function Representant(R,v,Phi::Vector,vals::Vector,Lambda::Vector)
  res=0
  A=base_ring(Phi[1])
  pi=BaseGenerators(A)
  k=length(pi)
  Gamma=map(valuation,pi)
  r=length(vals)
  F=parent(Lambda[1][1])
  n=degree(F) ## C'est quoi ca ? -> marche sur les corps finis mais pas pour tout visiblement. GENERIQUE
  maxdegs=[[degree(Phi[i])-1 for i in 2:r];degree(R)]
  for i in 0:degree(R)
    if coeff(R,i)!=0
      good_vals=AllMonomialsGivenV(maxdegs,vals,Gamma,v-i*vals[k])
      if length(good_vals)==1 res+=coeff(R,i)*prod([pi[j]^good_vals[1][j] for j in 1:k])*prod(Phi[i]^good_vals[i+1] for i in 1:r-1)
      else
        ## This code is not correct if A is defined over some extension field (like GF(9)) GENERIQUE
        # Building possible R-value for each "monomial" with the prescribed valuation.
        tmp=[prod([Lambda[1][j]^val[1][j] for j in 1:k])*prod([Lambda[j]^val[j] for j in 2:r]) for val in good_vals]
        # matrix passage between the basis generated by the R-value and the one of the residue field.
        M=matrix(ResidueField(A),n,n,[coeff(tmp[j],l) for l in 0:n-1 for j in 1:n])
        target=[coeff(coeff(R,i),j) for j in 0:n-1]
        co=inv(M)*target
        toto=sum([co[l]*prod([pi[j]^Int(good_vals[l][1][j]) for j in 1:k])*prod([Phi[j-1]^Int(good_vals[l][j]) for j in 2:r])*Phi[r]^i for l in 1:n])
        res+=toto
      end
    end
  end
  return res
end
