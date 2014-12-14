
function p = Zp(a, b)

if nargin == 0
   p.c = [];
   p.p = 30097;
   p = class(p,'Zp');

elseif isa(a,'Zp')

   p = a;

else

	if nargin == 2

		
		if isa(b,'Zp')
			prime = b.p;
		else
			prime = b;
		end
	else
		prime = 30097;
	end

   p.c = mod(a, prime);
   p.p = prime;
   p = class(p,'Zp');
end
  