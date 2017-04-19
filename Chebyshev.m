%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class Chebyshev                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To do:
%   o I would eventually like to incorporate the concept of a null entry
%     that could work naturally with times() and disp() -- possibly using 
%     Chebyshev(0,'A',0) -- to get rid of the isempty() checks.
%   o Enforce reservation of 'A' variable
%
% Changes:
%   9/7/2011 + changed all +'s into horzcats in diff() to avoid calling 
%              combine() every iteration
%            + added prune() and reduce()
%  9/20/2011 + corrected diff()
%            + inner() now works for sums and products; input can be one
%              or two arguments
%            + added a check in horzcat to omit zeros in sums
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Chebyshev < handle
    
    properties (SetAccess = private)
        coeff
        degree
        symbol
    end
    
    methods
        function T = Chebyshev(varargin)
           switch(size(varargin,2))
                % (degree, variable, coefficient)
                case 3  
                    T.degree = varargin{1};
                    s = varargin{2};

                    if T.degree ~= 0
                        disp('Can only set coefficient for T_0')
                        return;
                    end

                    T.coeff = varargin{3};

                    if ischar(s)
                        T.symbol = sym(s);
                    else
                        T.symbol = s;
                    end
                    
                % (degree, variable) or (list of coeffs, variable)
                case 2
                    
                    % More than one element in list -> list of coeffs
                    if length(varargin{1}) > 1
                        T = 0;
                        b = varargin{1};
                        c = length(varargin{1});
                        d = varargin{2};
                        
                        for i = 1:length(varargin{1})-1
                            T = T + b(i)*Chebyshev(c-1,d);
                            c = c - 1;
                        end
                        
                        T = T + b(i+1)*Chebyshev(c-1,'A');
                        
                    % Just one element -> degree
                    else
                        T.coeff = 1;
                        T.degree = varargin{1};
                        s = varargin{2};

                        if ischar(s)
                            T.symbol = sym(s);
                        else
                            T.symbol = s;
                        end
                    end
                
                % Another Chebyshev polynomial
                case 1
                    if ~isa(varargin{1},'Chebyshev')
                        disp('Argument must be an instance of the class.')
                        return;
                    end
                    
                    c = varargin{1};
                    T = c;
           end % switch
        end % constructor
        
        function disp(obj)
            disp(tag(obj))
        end % disp
        
        function r = inner(varargin)
        % Inner product
        
            if size(varargin,2) == 2
                
                a = varargin{1};
                b = varargin{2};
                
                if ~isa(a,'Chebyshev') || ~isa(b,'Chebyshev')
                    disp('Inner: both operands must be objects of the Chebyshev class.')
                    return
                end

                c = a*b;
                r = 0;

                % Iterate through terms in the sum
                for i = 1:size(c,2)
                    d = reduce(c(:,i));
                    r = r + inner(d);   
                end
                
            elseif size(varargin,2) == 1 && size(varargin{1},2) > 1
                
                r = 0;
                
                c = varargin{1};

                % Iterate through terms in the sum
                for i = 1:size(c,2)
                    d = reduce(c(:,i));
                    r = r + inner(d);   
                end
                
            elseif size(varargin,2) == 1 && size(varargin{1},2) == 1
                
                d = varargin{1};
                
                % Separate out terms with like variables
                M = zeros(length(d));
                M(1,1) = 1;
                k = 1;

                for j = 2:length(d)
                    if d(j).symbol ~= d(j-1).symbol 
                        k = k + 1;
                    end

                    M(j,k) = 1;
                    
                    if sum(M(:,k)) > 2
                        disp('Inner: too many terms of one variable.')
                        return
                    end
                end
                
                r = 1;
                
                for j = 1:size(M,2)
                    if sum(M(:,j)) == 0
                        break
                    end
                    
                    c = d.mask(M(:,j));
                    
                    if size(c,1) == 1 && c.symbol == 'A'
                        r = r*c.coeff;
                    elseif size(c,1) == 1 && c.symbol ~= 'A'
                        r = 0;
                        break
                    else
                    
                        % Assuming each product contains only two terms
                        if c(1).degree == c(2).degree
                            if c(1).degree == 0
                                r = r*pi;
                            else
                                r = r*pi/2;
                            end
                        else
                            r = 0;
                            break
                        end
                    end % if
                end
            elseif size(varargin,2) > 2
                disp('Inner: too many arguments.')
                return
            end % if
            
        end % inner

        function R = simplify(obj)
            
            % In the case that it's just a number, no need to simplify
            if size(obj,1) == 1 && size(obj,2) == 1
                R = obj;
                return
            end
            
            R = [];
            
            % Iterate through each term in the sum
            for i = 1:size(obj,2)
                
                % DEBUG
                %disp(['Simplifying ', tag(obj(:,i))]);
                
                q = reduce(obj(:,i));
                M = tree(q);
                s = Chebyshev(0,'A');
                
                % Iterate through each column of the tree matrix
                for j = 1:size(M,2)
                    m = M(:,j);
                    
                    % If all entries are 0, we're finished
                    if sum(m) == 0
                        break
                    end
                    
                    r = q.mask(m);
                    s = s*simp(r);
                end
                
                % Check whether this term requires further simplification
                recurse = 0;
                for j = 1:size(s,2)
                    n = sum(tree(s(:,j)));
                    
                    for k = 1:length(n)
                        if n(k) > 1
                            recurse = 1;
                        end
                    end
                end
                
                if recurse
                    s = simplify(s);
                end
                
                if isempty(R)
                    R = s;
                else
                    R = R + s;
                end
            end % for
           
        end % simplify
        
        function r = mtimes(a,b)
            r = [];
            
            if isa(a,'Chebyshev') && isa(b,'Chebyshev')

                for i = 1:size(a,2)
                    for j = 1:size(b,2)

                        if isempty(r)
                            r = [a(:,i);b(:,j)];
                        else
                            r = [r [a(:,i);b(:,j)]];
                        end

                    end
                end
                
                R = [];
                for i = 1:size(r,2)

                    % Put coefficients in front
                    r(:,i) = sort(r(:,i));

                    j = 1; 
                    q = Chebyshev(0,'A');

                    % Consolidate coefficients into one symbol
                    while j <= size(r,1) && r(j,i).symbol == 'A'
                        q = Chebyshev(0,'A',q.coeff*r(j,i).coeff);
                        j = j + 1;
                    end

                    % If this is the one and only term
                    if j == 1
                        if isempty(R)
                            R = r(:,i);
                        else
                            R = [R r(:,i)];
                        end

                    % Otherwise
                    else
                        if q.symbol == 'A' && q.coeff == 0
                            R = Chebyshev(0,'A',0);
                        else
                            if isempty(R)
                                R = [q; r(j:end,i)];
                            else
                                R = [R [q; r(j:end,i)]];
                            end
                        end
                    end
                end % for
                
                if isempty(R)
                    r = Chebyshev(0,'A',0);
                else
                    r = R;
                end
                
            elseif (isa(a,'Chebyshev') && isnumeric(b))
                
                if b == 0
                    r = Chebyshev(0,'A',0);
                elseif b == 1
                    r = a;
                else
                    r = Chebyshev(0,'A',b)*a;
                end
                
            elseif (isa(b,'Chebyshev') && isnumeric(a))
                
                if a == 0
                    % r = 0;
                    r = Chebyshev(0,'A',0);
                elseif a == 1
                    r = b;
                else
                    r = Chebyshev(0,'A',a)*b;
                end
                
            end
        end % times
        
        function r = plus(a,b)
        % Overloaded addition operator
        
            if isa(a,'Chebyshev') && isa(b,'Chebyshev')      
                r = [a b];
            elseif (isa(a,'Chebyshev') && isnumeric(b))
                if b == 0
                    r = a;
                else
                    r = plus(Chebyshev(0,a.symbol,b), a);
                end
            elseif (isa(b,'Chebyshev') && isnumeric(a))
                if a == 0
                    r = b;
                else
                    r = plus(Chebyshev(0,b.symbol,a), b);
                end
            end
            
            if isempty(r)
                r = Chebyshev(0,'A',0);
            else
                r = combine(r);
            end
        end % plus
        
        function r = prune(obj,prec)
        % Removes terms whose coefficients are smaller than prec
        
            r = Chebyshev(0,'A',0);
            for i = 1:size(obj,2)
                if obj(1,i).symbol == 'A' && abs(obj(1,i).coeff) > prec
                    r = r + obj(:,i);
                end
            end
        end % prune
        
        function R = combine(obj)
        % Combines like terms
            Q = ones(size(obj,2),1);
            k = 1;
        
            % Column
            for i = 1:size(obj,2)
                % Row
                for j = 1:size(obj,1)
                    
                    % Compute hash
                    Q(i) = Q(i)*bitxor(double(char(obj(j,i).symbol)),obj(j,i).degree);
                    
                end
                k = k + 1;
            end
            
            k = 2;
            n = 1;
            R = [];
            done = zeros(1,1);
            
            % Iterate through columns
            for i = 1:length(Q)
                
                % Check that the term hasn't been picked up already
                for l = 1:length(done)
                    if done(l) == i
                        B = 0;
                        break;
                    else
                        B = 1;
                    end
                end
                
                % Not picked up
                if B
                    
                    D = obj(:,i);
                    
                    for j = k:length(Q)

                        % Found a match
                        if (Q(i) == Q(j))
                            
                            % Coefficient of pivot term
                            if D(1).symbol == 'A'
                                c = D(1).coeff;
                            else
                                c = 1;
                            end
                            
                            % Coefficient of current term
                            if obj(1,j).symbol == 'A'
                                c = c + obj(1,j).coeff;
                            else
                                c = c + 1;
                            end
                            
                            % Strip off coefficient and replace
                            if D(1).symbol == 'A'
                                D = D(2:end);
                            end
                            
                            D = c*D;

                            % Log it
                            done(n) = j;
                            n = n + 1;
                            
                        end % if
                    end % for
                    
                    R = [R D];
                
                end % if
                k = k + 1;
            end % for
            
        end % combine
        
        function r = sympoly(obj)
            R = sym(zeros(size(obj)));
            
            for i = 1:size(obj,1)
                for j = 1:size(obj,2)
                    if char(obj(i,j).symbol) == 'A'
                        R(i,j) = obj(i,j).coeff;
                    else
                        R(i,j) = poly2sym(obj(i,j).poly,obj(i,j).symbol);
                    end
                end
            end
            
            r = 0;
            for i = 1:size(R,2)
                s = 1;
                for j = 1:size(R,1)
                    s = s*R(j,i);
                end
                r = r + s;
            end
            
            r = expand(r);
        end % sympoly
        
        function r = horzcat(a,b)
            
            % Check whether argument is null
            if size(a,1) == 0 || size(a,2) == 0 || (size(a,1) == 1 && a(1).symbol == 'A' && a(1).coeff == 0)
                r = b;
            elseif size(b,1) == 0 || size(b,2) == 0 || (size(b,1) == 1 && b(1).symbol == 'A' && b(1).coeff == 0)
                r = a;
            else
                
                % Check whether one column vector is longer than the other
                if size(a,1) > size(b,1)
                    
                    k = size(b,1);
                    for i = 1:size(a,1)-size(b,1)
                        for j = 1:size(b,2)
                            
                            % Pad the vector to make it fit
                            b(k+i,j) = Chebyshev(0,'A');
                        end
                    end
                elseif size(b,1) > size(a,1)
                    
                    k = size(a,1);
                    for i = 1:size(b,1)-size(a,1)
                        for j = 1:size(a,2)
                            
                            % Pad the vector to make it fit
                            a(k+i,j) = Chebyshev(0,'A');
                        end
                    end
                end
                
                % Reform the matrix
                r = Chebyshev(0,'A');
                for i = 1:size(a,1)
                    for j = 1:size(a,2)
                        r(i,j) = a(i,j);
                    end

                    for j = 1:size(b,2)
                        r(i,j+size(a,2)) = b(i,j);
                    end
                end
                
            end % if
        end % horzcat
        
        function r = diff(obj,var)
            r = [];
            
            % Iterate through each term in a sum
            for k = 1:size(obj,2)
                
                % Iterate through each term in a product
                for i = 1:size(obj,1)
                    s = [];

                    % Is the term a function of differentiation variable?
                    if char(obj(i,k).symbol) == var
                        n = obj(i,k).degree;
                        sym = obj(i,k).symbol;

                        % Form derivative
                        for j = 0:n-1
                            if mod(n,2)
                            % n even
                                if isempty(s)
                                    s = Chebyshev(2*j+1,sym);
                                else
                                    s = [s Chebyshev(2*j+1,sym)];
                                end
                            else
                            % n odd
                                if isempty(s)
                                    s = Chebyshev(2*j,sym);
                                else
                                    s = [s Chebyshev(2*j,sym)];
                                end
                            end
                        end
                        
                        if mod(n,2)
                            s = s + Chebyshev(0,'A',-1);
                        end

                        if size(obj,1) > 1
                            % First term in product
                            if i == 1
                                if isempty(r)
                                    r = 2*n*s*obj(i+1:end,k);
                                else
                                    r = [r 2*n*s*obj(i+1:end,k)];
                                end

                            % Last term in product
                            elseif i == size(obj,1)
                                if isempty(r)
                                    r = 2*n*s*obj(1:i-1,k);
                                else
                                    r = [r 2*n*s*obj(1:i-1,k)];
                                end

                            % Middle terms in product   
                            else
                                if isempty(r)
                                    r = 2*n*s*obj(1:i-1,k)*obj(i+1:end,k);
                                else
                                    r = [r 2*n*s*obj(1:i-1,k)*obj(i+1:end,k)];
                                end
                            end
                        else
                            if isempty(r)
                                r = 2*n*s;
                            else
                                % r = r + 2*n*s;
                                r = [r 2*n*s];
                            end % if
                        end % if
                    end % if
                end % for
            end % for
            if isempty(r)
                r = Chebyshev(0,'A',0);
            else
                % Combine like terms
                r = combine(r);
            end
        end % diff
        
        function [obj,idx,varargout] = sort(obj,varargin)
        % Overloaded sort method    
        
            varargout = cell(1,nargout-2);
            [~,idx,varargout{:}] = sort([obj.symbol],varargin{:});
            obj=obj(idx);
        end % sort
        
        function r = reduce(obj)
        % Removes extra rows for further processing
        
            if size(obj,2) > 1
                disp(['reduce(): can only reduce products. (', tag(obj), ')'])
                r = [];
            else
                c = 1;
                for i = 1:length(obj)
                    if obj(i).symbol ~= 'A'
                        break
                    else
                        c = c*obj(i).coeff;
                    end
                end

                r = c*obj(i:end);
            end
                
        end % reduce
    
    end % public methods
    
    methods (Access = private)
        
        function M = tree(obj)
            s = size(obj,1);
            M = zeros(s,s+1);
            
            obj = sort(obj);
            
            k = 1;
            q = obj(1,1).symbol;
            M(1,1) = 1;
            
            for i = 2:size(obj,1)
                
                if obj(i,1).symbol == q && sum(M(:,k)) ~= 2
                    M(i,k) = 1;
                elseif obj(i,1).symbol ~= q || sum(M(:,k)) == 2
                    k = k + 1;
                    M(i,k) = 1;
                end
                
                q = obj(i,1).symbol;     
            end
        end % tree
        
        function q = mask(obj,mask)
        % Input is a vector consisting of 1's or 0's, indicating which
        % terms to keep or zero, respectively
            
            q = [];
            obj = sort(obj);
            
            for i = 1:size(obj,2)
                 r = [];
                 
                for j = 1:size(obj,1)
                    if mask(j)
                        r = [r; obj(j,i)];
                    end
                end
                
                q = [q r];
            end
        end % mask
        
        function r = poly(obj)
            r = ChebT(obj.degree);
        end % poly
        
        function f = ChebT(n)
            if n == 0
                f = 1;
            elseif n == 1
                f = [1 0];
            else
                f = [2*ChebT(n-1) 0] - [0 0 ChebT(n-2)];
            end
        end % ChebT
        
        function r = tag(obj)
            
            % Product or sum of terms
            if size(obj,1) > 1 || size(obj,2) > 1
                r = [];
                
                % Iterate through columns
                for i = 1:size(obj,2)
                    
                    % If sum
                    if i > 1
                        r = [r ' + '];
                    end
                    
                    term = reduce(obj(:,i));
                    
                    for j = 1:length(term)

                        if j > 1
                            r = [r '*'];
                        else
                            r = r;
                        end

                        r = [r tag(term(j))];

                    end % for
                end
                
            % Single term
            else
                if obj.degree == 0
                    r = num2str(obj.coeff);
                else
                    r = ['T_' num2str(obj.degree) '(' char(obj.symbol) ')'];
                end
            end
        end % tag
        
        function r = simp(obj)
            if size(obj,1) < 2 || obj(1).symbol ~= obj(2).symbol
                r = obj;
            else
                a = obj(1).degree;
                b = obj(2).degree;
                s = obj(1).symbol;

                r = 1/2*Chebyshev(a+b,s) + 1/2*Chebyshev(abs(a-b),s);
            end
        end % simp
            
    end % private methods
    
end % classdef