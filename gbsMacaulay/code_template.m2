-- Macaulay2 code template for gbsolver generator
-- 
-- by Martin Bujnak, sep2008
--

KK = ZZ/$PRIME$

R = KK[x_1..x_$UNKCNT$, MonomialOrder=>GRevLex]

-- generated part
$EQUATIONS$
$EQUATIONSVEC$
-- generated part end

I1 = ideal(f); 
gbTrace 3;
dm = dim I1;
dg = degree I1;

-- Do not modify following lines
print "dim:"
print dm;
print "deg:" 
print dg;

A = R/I1;
Ab = basis A;

print "@@@@@@@@@@@@@@"
print Ab;
print "@@@@@@@@@@@@@@"

exit 0
