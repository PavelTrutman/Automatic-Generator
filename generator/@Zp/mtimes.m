
function r = mtimes(p,q)

	% Zp/PLUS  Implement (p + q) mod Zp.

	p = Zp(p, q);
	q = Zp(q, p);

	r = Zp(0, p);

	r.c = mod((p.c * q.c), p.p);
end
