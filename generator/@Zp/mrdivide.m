function [c] = mrdivide(p,q)

    p = Zp(p, q);
    q = Zp(q, p);

    % invert q
    [g,c,d] = gcd(q.c,-q.p);
    c = mod(c, q.p) ;

	c  = Zp(p.c*c, p.p);

end % mrdivide-method
