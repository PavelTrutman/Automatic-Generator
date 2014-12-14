
function r = mpower(p,q)

	% ZP/MPOWER   Implement p ^ q for Zp numbers.

	if q < 0

        [g,c,d] = gcd(p.c, -p.p);
        p.c = mod(c, p.p);
		q = -q;
	end

	r = Zp(1, p.p);
	for i=1:q
	 	
		r.c = r.c*p.c;
	end	
	r.c = mod(r.c, r.p);
end
