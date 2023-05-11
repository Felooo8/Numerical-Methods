function [c, fc, i] = bisection(f, a, b, delta)
%
%   BISECTION
%       Finds the root of a function of one variable
%       using bisection method
%  
%   INPUT
%       f      -  given function  
%       a      -  starting point
%       b      -  ending point
%       delta  -  accuracy
%
%   OUTPUT
%       c     -  root 
%       fc     -  function value in point c f(c)
%       i   -  number of iterations
%
%   SAMPLE RUN
%       >> [c, fc, i] = bisection(@(x) 2.2*x.*cos(x)-2*log(x+2), 2, 7, 1e-8)
%       ans =
%           5.0707
%
    
    if nargin < 4
        delta = 1e-8;
    end

    fc=inf; % to run while loop
    i=0;
    while (abs(fc) > delta && abs(b-a)>1e-8)
        i=i+1;

        c=(a+b)/2; % gets points c between a and b
        fc=f(c);
        fa=f(a);

        if (fa*fc<0) % if root is between a and c
            b=c;
        else
            a=c;
        end
    end
end