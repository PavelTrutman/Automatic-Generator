% reorder a set of monomials stored in cell array.
%
% function [mons] = ReorderMonomials(P, unknowns, orderingfnc)
%
%   P - set of monomials
%
%   unknowns are list of unknown variables (reminding vars are considered
%   as const) - example unknowns = {'x_0' 'x_1' 'x_2' 'x_3' 'x_4' 'x_5')
%
%   orderingfnc(a,b) compare two monomials (where a, b are sets of powers)
%   and return 1 to indicat a>b or -1 to indicat a<b. 
%   orderingfnc is pointer to a function, however few presets are available
%   set:  
%        1 - lex ordering
%        2 - graded lex ordering
%        3 - graded reverse lex ordering
%   
%
% warning! unknowns are of the form 'x_0' ... 'x_N' or completely different
%          letters (not mixtures! like 'a' 'x_4' 'b' ... )% 
%
% by Martin Bujnak, nov2007
% last edit by Pavel Trutman, oct 2014

function [monomials degree] = ReorderMonomials(P, unknowns, orderingfnc)

    monomials = [];
    degree = [];
    unkcnt = length(unknowns);

    % fix degrees and local ordering + remove duplicates
    for i = 1:length(P)
        p = P{i};
        [monomial mcnt]= Monomial_PerformSubst(unknowns, {}, p);
        monomial = Monomials_order(monomial, unkcnt+1, 'x_1');
        monomials = union(monomials, cellstr(monomial));
    end      
    
   % sort monomials globaly
    monscnt = size(monomials, 2);
    monsdeg = zeros(1, monscnt);
    monsdegre = zeros(monscnt, unkcnt);
    for i=1:monscnt 

        [monsdeg(i) monsdegre(i, :)]= Monomial_degree(monomials(i), unkcnt);
     
        % back subsitute monomials...
        monomials(i) = Monomial_PerformSubst({}, unknowns, monomials(i));
    end
    
    if nargin > 2
        
        if isnumeric(orderingfnc)
            
            switch orderingfnc
                case 1
                    orderingfnc = @Monomial_lexorder;
                case 2
                    orderingfnc = @Monomial_grlexorder;
                case 3
                    orderingfnc = @Monomial_grevlexorder;
            end
        end
        
        % reorder monomials (max-sort)
        for i=1:size(monomials, 2)
            
            maxmon = monsdegre(i,:);
            maxidx = i;
            for j=i+1:size(monomials, 2)
                
                if orderingfnc(monsdegre(j,:), maxmon) > 0
                    
                    maxmon = monsdegre(j,:);
                    maxidx = j;
                end
            end
            
            % swap monomials
            zl = monomials(i);
            monomials(i) = monomials(maxidx);
            monomials(maxidx) = zl;
            
            monsdegre(maxidx,:) = monsdegre(i,:);
            monsdegre(i,:) = maxmon;
            
            degree(i,:) = maxmon;
        end
        
    else
        [v i] = sort(monsdeg,'descend');
        monomials = monomials(i);
        degree = monsdegre(i, :);
    end    
end