% monomials extraction utility 
% (c) Martin Bujnak, martinb@dataexpert.sk, aug2006
%
% [monomials degree] = ExtractMonomials(P, unknowns, orderingfnc)
% usage:
%   P is cell array of polynomial equations
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
%          letters (not mixtures! like 'a' 'x_4' 'b' ... )
%
% 
function [monomials degree] = ExtractMonomials(P, unknowns, orderingfnc)

    uns = '#';
    for i = 1:size(unknowns, 2)

        uns = [uns char(unknowns{i}) '#'];
    end

    unkcnt = size(unknowns, 2);
    monomials = [];
    degree = [];
    w = size(P, 2);
    for i = 1:w

        if ~iscell(P) 
            
            Pi = expand(P(i));
        else
            
            Pi = expand(P{i});
        end

        % extract elems
        elems = strread(char(Pi), '%s',  'delimiter', '+-');

        % extract unknowns
        for j = 1:size(elems, 1)

            if ~isempty(char(elems(j)))

                monomial = '';
                tvars = strread(char(elems(j)), '%s',  'delimiter', '*');
                divs = [];
                vars = [];

                for k = 1:size(tvars,1)

                    dvars = strread(char(tvars(k)), '%s',  'delimiter', '/');
                    vars = [vars; dvars(1)];
                    
                    if size(dvars,2) > 1
                        divs = [divs; dvars(2:end)];
                    end
                end
                
                div = size(vars, 1);
                if ~isempty(divs)
                    vars = [vars; divs];
                end
                
                for k = 1:size(vars,1)

                    e = strread(char(vars(k)), '%s',  'delimiter', '^');
                    fs = ['#' strtrim(char(e(1))) '#'];

                    if ~isempty(strfind(uns, fs))

                        if isempty(monomial)
                            monomial = [monomial strtrim(char(vars(k))) ];
                        else

                            if k > div
                                monomial = [monomial '/' strtrim(char(vars(k))) ];
                            else
                                monomial = [monomial '*' strtrim(char(vars(k))) ];
                            end
                        end
                    end

                end

                if ~isempty(monomial)

                    % Monomials_order requires monomial singletons to be of
                    % the form "x_1 ... x_N" - this defines ordering too.
                    % First we substitute monomials and then fix it back.
                    
                    [monomial mcnt]= Monomial_PerformSubst(unknowns, {}, monomial);
                    monomial = Monomials_order(monomial, unkcnt+1, 'x_1');

                    monomials = union(monomials, cellstr(monomial));
                end
            end
        end

    end

    % sort monimials by degree (+ lex ordering)
    monscnt = size(monomials, 2);
    monsdeg = zeros(1, monscnt);
    monsdegre = zeros(monscnt, unkcnt);
    for i=1:monscnt 

        [monsdeg(i) mdeg]= Monomial_degree(monomials(i), unkcnt);
        monsdegre(i, :) = mdeg;
     
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
                
                if orderingfnc(maxmon, monsdegre(j,:)) < 0
                    
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
