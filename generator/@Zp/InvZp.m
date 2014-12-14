% Zuzana Kukelova 2007-08-21
% Inverse element to the element from Zp
% Input:    x - element from Zp
%           p - prime
%           
% Output: c = x^(-1) in Zp


function c =  InvZp(x,p)

% return x*c - d*p = g 
% since p is prime number g=1 => c = x^(-1)
[g,c,d] = gcd(x,-p);

c = mod(c,p);
