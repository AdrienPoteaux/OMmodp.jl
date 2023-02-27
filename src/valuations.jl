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
  # comparison function
  function cmp(A,B)
    O = L[1]
    turn = (A[1]-O[1])*(B[2]-O[2]) - (A[2]-O[2])*(B[1]-O[1])
    if turn == 0 # A, 0 ,B on the same line.
      return A[1] > B[1]
    end
    return turn > 0 # turn > 0 means OA < OB
  end
  sorted = [[L[1]] ; sort(L[2:length(L)],lt=cmp)]
  s=Stack{Vector{Rational}}()
  push!(s,L[1])
  push!(s,L[2])
  i=3
  while i<length(L)
    if length(s) == 1
      push!(s,L[i])
      i+=1
      continue
    end
    curr = first(s)
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

"""
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
  return minimum(map(valuation,F.coeffs))
end

"""
    PhiVal(elt::Vector, vals::Vector)

Compute the valuation of elt

# arguments
- elt a table representing some Phi-adic expansion (i.e. a table of table of ... of polynomials), as given by PhiExp
- vals the (augmented) valuations of the Phi[i]
"""
function PhiVal(elt::Vector, vals::Vector)
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
  #=
  Remark : we might want to change the implementation to get the "classical" Newton polygon
  =#
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
  # Now storing the points we have to consider for the lowest convex hull
  left_points = [[first-1,valuations[first]]] # corresponds to negative slopes
  if first != l
    tmp = (valuations[l]-valuations[first])//(l-first)
    a = -numerator(tmp)
    b = denominator(tmp)
    c = a*(first-1) + b*valuations[first]
    # the line from [first-1,valuations[first]] to [l-1,valuations[l]] is a*i+b*j=c
    for i in first+1:l-1
      if (a*(i-1)+b*valuations[i]<c) # this point can be in the lowest convex hull.
        left_points = [left_points;[[i-1,valuations[i]]]]
      end
    end
    left_points = [left_points;[[l-1,valuations[l]]]]
  end
  points = LowerConvexHull(left_points)
  if l!=r
    points = [points;[[r-1,valuations[r]]]]
  end
  if r != last
    right_points=[[r-1,valuations[r]]]
    tmp = (valuations[r]-valuations[last])//(r-last)
    a = -numerator(tmp)
    b = denominator(tmp)
    c = a*(r-1) + b*valuations[r]
    # the line from [r-1,valuations[r]] to [last-1,valuations[last]] is a*i+b*j=c
    for i in r+1:last-1
      if (a*(i-1)+b*valuations[i]<c) # this point can be in the lowest convex hull.
        right_points=[right_points;[[i-1,valuations[i]]]]
      end
    end
    right_points=[right_points;[[last-1,valuations[last]]]]
    tmp = LowerConvexHull(right_points)
    points=[points;tmp[2:length(tmp)]]
  end
  return points
end

### ICI ICI ICI
# More an intern function than anything else
# GENERIQUE : il a pas l'air d'aimer le elt::Vector (je l'ai vire du coup)
function AllCoeffGivenV(elt, vals::Vector, v)
#  In: elt a Phi expansion of some polynomial, vals the valuation of the key polynomials and v some targeted valuation
# Out: all coefficients having valuation v, taking into account all phi
  k = length(vals)
  res=[]
  if k>0
    for i in eachindex(elt)
      tmp = AllCoeffGivenV(elt[i], vals[1:k-1], v-(i-1)*vals[k])
      for j in eachindex(tmp)
        tmp[j] = [tmp[j] ; i-1]
      end
      if length(tmp)>0#necessaire ?
        res = [res ; tmp]
      end
    end
    return res
  end
  # k == 0
  res=[]
  if denominator(v)!=1
    return res
  end
  for i in 1:length(elt)
    res=[res;[[coeff(coeff(elt,i-1),numerator(v)),numerator(v),i-1]]] # NON GENERIQUE ?
  end
  return res
end

#=
On the formulas to compute the residual polynomial

We have R_0 = (F/x^v_0(F))(0,z)

Issac version : R_k(F) = \sum_i z_{k-1}^\tau_{k,i}*R_{k-1}(a_i*phi_k^i)(z_{k-1}) * z^{(i-i_k(F))/q_k}

\tau_{k,i} = (i_{k-1}(a_i)+\beta_{k-1}*v_k(a_i*phi_k^i))/q_{k-1}

One example :

P = y^2 - 2*x

On veut (T^2+T+1)^2 = T^4 + 2*T^3 + 3*T^2 + 2*T +1 comme polynome caracteristique

F = P^4 + 2*x*y*P^3 + 6*x^3*P^2 + 4*x^4*y*P + 4*x^6

Puiseux fait x=2*x^2, y=x*(y+2) ; donc P devient 4*y+...

Du coup le polynome de facette devient 2^8*y^4 + 2^9*y^3 + 3*2^8*y^2 + 2^9*y + 2^8, ce qu'on veut.

polynome residuel (formules Montes) : beta_1 = 1 ; v_2(x)=2, v_2(y)=1, v_2(P)=2
a_0 = 4*x^6   -> R_1(a_0)=4 ; \tau_{1,0} = (0+1*v_2(4*x^6))/2 = 6     --> 2^8
a_1 = 4*x^4*y -> R_1(a_1)=4 ; \tau_{1,1} = (1+1*v_2(4*x^4*y*P))/2 = 6 --> 2^8*z
a_2 = 6*x^3   -> R_1(a_2)=6 ; \tau_{1,2} = (0+1*v_2(6*x^3*P^2)) = 5   --> 3*2^6*z^2
a_3 = 2*x*y   -> R_1(a_3)=2 ; \tau_{1,3} = (1+1*v_2(2*x*y*P^3))/2 = 5 --> 2^6*z^3
a_4 = 1       -> R_1(a_4)=1 ; \tau_{1,4} = (0+v_2(P^4))/2 = 4         --> 2^4*z^4

Je trouve donc 2^4*(z^4+4*z^3+12*z^2+16*z+16)=2^4*(z^2+2*z+4)^2

R_2(2*z) est egal au polynome caracteristique !

Quid avec des valuations "rationnelles" ?
To deal with:
  1) on a plus vraiment de m et q, donc de beta.
  2) les valeurs des valuations sont differentes.


=#
"""
    PhiResidualPol(elt::Vector, Delta, vals::Vector, Lambda::Vector)

Compute the residual polynomial of elt according to the edge Delta.

vals provide the list of augmented valuations
Lambda the list of "correcting terms"
"""
function PhiResidualPol(elt::Vector, Delta, vals::Vector, Lambda::Vector)
  slope = (Delta[2][2]-Delta[1][2])/(Delta[2][1]-Delta[1][1])
  m=-numerator(slope)
  q=denominator(slope)

  for i in Delta[1][1]:q:Delta[2][1]
    j = numerator((i-Delta[1][1])//q) # q divides i-Delta[1][1]
    tmp = AllCoeffGivenV(elt[i],vals,Delta[1][2]-(j*m))

  end
end
