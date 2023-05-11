function [R] = laguerre_with_horner(a,X,delta)
%
% Function that finds roots of given polynomial using Laguerre and Horner
% methods (defined below)
%  
%   INPUT
%       a     -  coefficients of a polynomial (desceding by degree)
%       X     -  vector of x0 for each iteration
%       delta  -  accuracy
%
%   OUTPUT
%       R     -  vector of roots 
%
%   SAMPLE RUN
%       >> [R] = task2([-2,3,3,1,2],[-1,2,2,2],1e-8)
%       ans =
%
%           2.3197 + 0.0000i
%           0.0419 - 0.6895i
%           0.0419 + 0.6895i
%           -0.9035 + 0.0000i
%

    if nargin < 2
        X = [1,1,1,1];
        delta = 1e-8;
    elseif nargin < 3
        delta = 1e-8;
    end

    n = size(a,2)-1;
    R=zeros(n,1); % initialize vector of roots
    iters=zeros(n,3); % matrix of: number of iterations, root, f(root)

    i=1;
    while i < n+1
        [root, y, iter] = laguerre(a,X(i),delta);

        R(i) = root;
        iters(i,:) = [iter,root,y];

        a = horner(a,root);
        
        % looks for conjugate of complex root
        im = imag(root);
        if im~=0
            conj_root = conj(root);

            a = horner(a,conj_root);

            i = i+1;
            R(i) = conj_root;
            iters(i,:) = [iter,conj_root,y];
        end

        i=i+1;
    end
end


function [x, y, i] = laguerre(f, x0, delta)
%
%   LAGUERRE
%       Finds the root of a function of one variable
%       using laguerre method
%  
%   INPUT
%       f      -  coefficients of a polynomial f(x) (desceding by degree)
%       x0     -  value of x for the initial iteration
%       delta  -  accuracy
%
%   OUTPUT
%       x     -  root 
%       y     -  function value in point c - f(c)
%       i     -  number of iterations
%
%   SAMPLE RUN
%       >> [x] = laguerre([-2,3,3,1,2], 2, 1e-8)
%       ans =
%           [2.3197,1.061335463958812e-10,3]
%

    n = size(f,2)-1;
    df = polyder(f);  %first differitial of f
    d2f = polyder(df); %second differitial of f
    X(1)=x0;
    i=1;
    
    while (abs(polyval(f,X(i))) > delta)
        yf=polyval(f,X(i));
        ydf=polyval(df,X(i));
        yd2f=polyval(d2f,X(i));
    
        G=sqrt((n-1)*((n-1)*ydf^2-n*yf*yd2f));
        
        D1=ydf-G;
        D2=ydf+G;
        if abs(D1) > abs(D2) % chooses bigger absolute value
            D=D1;
        else
            D=D2;
        end
    
        X(i+1)=X(i)-n*yf/D; %#ok<AGROW> 
        i=i+1;
    end
    x=X(i);
    y=polyval(f,x);
end


function [Q] = horner(a,root)
%
%   HORNER
%       Returns coefficients of Q(x), where f(x)=(x-root)Q(x)
%       using Horner's method
%  
%   INPUT
%       a      -  coefficients of a polynomial f(x) (desceding by degree)
%       root   -  root of f(x) by which we divide
%
%   OUTPUT
%       q      -  coefficients of a new polynomial, where f(x)=(x-root)Q(x)
%
%   SAMPLE RUN
%       >> [q] = horner([1,2,1], -1)
%       ans =
%           [1,1]
%
    
    n = length(a)-1; % degree of Q(x)
    q=zeros(1,n-1);
    q(1) = a(1);
    for i = 2:n
        q(i) =  a(i) + q(i-1)*root;
    end
    Q=q;
end