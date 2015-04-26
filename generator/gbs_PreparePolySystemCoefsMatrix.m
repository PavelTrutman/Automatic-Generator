% Generate all polynomials required to build an action matrix for
% variable "actMvar" (detect actMvar if variable was not specified)
% (GBsolver subroutine)
% by Martin Bujnak, mar2008
% last edit by Pavel Trutman, March 2015


function [M, trace, symcoefs, amVar, amLT, amLTall, algBidx, algB] = gbs_PreparePolySystemCoefsMatrix(cfg, eq, known, unknown, algB, amVar)

    ordering = cfg.ordering;

    % create random instance
    if ~isfield(cfg, 'eqinstance')

        cfg.eqinstance = cfg.InstanceGenerator(cfg, eq, known, unknown);
    end
    
    if nargin < 5 || isempty(algB)
        
        % find quotion ring basis ( monomial basis of quotient ring C[x1,...xn]/I )
        [algB] = cfg.GBSolver(cfg, eq, known, unknown);
        
        if isempty(algB)
            
            M=[];
            trace=[];
            symcoefs=[];
            amVar=[];
            amLT=[];
            amLTall=[];
            algBidx=[];
            return;
        end
    end
    
    % prepare algebra B and find which monomials we need to build the
    % action matrix
    fprintf('analyzing quotient ring basis\n');
    
    if nargin < 9
        
        amVars = unknown;
    else
        amVars = amVar;
    end
    
    amStats = {};
    for i=1:length(amVars)
        
        actMvar = amVars{i};
        
        [amLT, amLTcnt, amLTall, algB, algBidx] = gbs_ParseAlgebraB(algB, actMvar, unknown);
        
        amStats{i}.actMvar = actMvar;
        amStats{i}.amLT = amLT;
        amStats{i}.amLTcnt = amLTcnt;
        amStats{i}.amLTall = amLTall;
        amStats{i}.algB = algB;
        amStats{i}.algBidx = algBidx;
    end

    %
    % extract equations into coeffients + monomial form
    fprintf('extracting coefficient & monomials\n');
    [monomials, p, ~, symcoefs, maxdeg] = gbs_ParseEquations(cfg, eq, known, unknown);
    
    
    %
    % just statistics (can be removed from this code)
    [monomials] = ReorderMonomials(monomials, unknown, ordering);
    
    fprintf('...used %d monomials : \n   ', length(monomials));
    for mons = monomials
        fprintf('%s ', char(mons));
    end
    fprintf('\n');

    fprintf('...max. poly. degree %d\n', maxdeg);

    
    if ~isempty(cfg.crashlog)

        fprintf('...crash log file detected - recovering from previous state using: %s\n', cfg.crashlog);
        [~, ~, max_deg] = gbs_ResetFromLog(cfg.crashlog);

        maxdeg = max_deg;
        fprintf('...need to generate polynomials to degree %d\n', maxdeg);
    end
    

    
    %
    % generate monomials up-to maxdeg (matrix cols)
    allmons = {'1'};
    alldegs = zeros(1, length(unknown));
    allmonsdeg = 0;
    for i=1:maxdeg
        [mons, degs] = GenerateMonomials(i, unknown, ordering);
        allmons = [mons allmons];
        alldegs = [degs; alldegs];
        allmonsdeg = [i*ones(1, length(mons)) allmonsdeg];
    end
    
    % create action var. masks
    cols = size(allmonsdeg, 2);
    for a=1:length(amStats)
        amStats{a}.zero_el = setdiff(1:cols, (cols+1)-amStats{a}.algBidx);
    end
    
    if cfg.exportEqs
        
        fprintf('...exporting parsed equations');
        
        % debug - create coefficient matrix with all initial polynomials
        cols = size(allmonsdeg, 2);
        Minit = zeros(length(eq), cols);
        for i=1:length(eq)
            for j=1:p{i}.monscnt
                pmon2 = p{i}.deg(j,:);
                order = GetMonomialOrder(pmon2, unknown); 
                Minit(i, cols-(order-1)) =  p{i}.coefsIDX(j);
            end        
        end
        
        nonzero = find(sum(Minit) ~= 0);
        Mmons = Minit(:, nonzero);
        
        gbs_DBGExport(known, unknown, eq, symcoefs, allmons, Minit, Mmons);
    end

    fprintf('Adding polynomials\n');
    
    PolynomialsGenerator = str2func(['gbs_GeneratePolynomials_', cfg.PolynomialsGenerator]);
    [foundVar, M, trace] = PolynomialsGenerator(p, eq, unknown, maxdeg, alldegs, allmonsdeg, allmons, amStats, cfg, cfg.PolynomialsGeneratorCfg);
    
    if ~foundVar
        
        M = [];
        trace = [];
        
        amVar = [];
        amLT = [];
        amLTall = [];
        algBidx = [];

    else
        
        amVar = amStats{foundVar}.actMvar;
        amLT = amStats{foundVar}.amLT;
        amLTall = amStats{foundVar}.amLTall;
        algBidx = amStats{foundVar}.algBidx;
    end
end

%{
function [M, Mcoefs] = gbs_ReallocMatrix(M, Mcoefs, torows, tocols)

    crows = size(M, 1);
    ccols = size(M, 2);
    
    if crows ~= torows
        
        M = [M; zeros(torows - crows, ccols)];
        Mcoefs = [Mcoefs; zeros(torows - crows, ccols)];
    end
    
    if ccols ~= tocols
               
        % shift coefs
        newcols = zeros(torows, tocols - ccols);
        M = [newcols M];
        Mcoefs = [newcols Mcoefs];
    end
end
%}

%{
function [c] = InvZp(x, p)

    [g,c] = gcd(x,-p);
    c = mod(c,p);
end
%}
