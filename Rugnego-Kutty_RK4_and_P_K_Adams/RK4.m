function [y] = RK4(dx1,dx2,x1,x2,h,a,b)
% Returns y, which represents a list of points (x1,x2) calculated using
% 4th-order Runge-Jutta method based on provided differential equations
% (dx1,dx2)
%  
%   INPUT
%       dx1     -  first differential equation
%       dx2     -  seconds differential equation
%       x1      -  initial x1 (dx1(0))
%       x2      -  initial x2 (dx2(0))
%       h       -  step lenght
%       a       -  starting point
%       b       -  stopping point
%
%   OUTPUT
%       y       -  vector of points (x1,x2) 
%
%SAMPLE RUN 
%       >> [y] = RK4(@(x1,x2)x2 + x1*(0.5-(x1)^2-(x2)^2),@(x1,x2)-x1 + x2*(0.5-(x1)^2-(x2)^2),0.03,0.001,0.01,0,20)
%
%
    tic;
    t=a:h:b;
    halfstep=h/2;
    X1(:,1) = x1;
    X2(:,1) = x2;
    for i=1:(length(t)-1)
        k11=dx1(x1,x2);
        k21=dx2(x1,x2);
            
        k12=dx1(x1+halfstep*k11,x2+halfstep*k21);
        k22=dx2(x1+halfstep*k11,x2+halfstep*k21);

        k13=dx1(x1+halfstep*k12,x2+halfstep*k22);
        k23=dx2(x1+halfstep*k12,x2+halfstep*k22);

        k14=dx1(x1+h*k13,x2+h*k23);
        k24=dx2(x1+h*k13,x2+h*k23);

        x1=x1+(h/6)*(k11+2*(k12+k13)+k14);
        x2=x2+(h/6)*(k21+2*(k22+k23)+k24);

        X1(:,i+1)=x1;
        X2(:,i+1)=x2;
    end
    toc;
    y = [X1; X2];
    plot(X1,X2);
    %plot(t(1:n),X1(1:n));
    %plot(t(1:n),X2(1:n));
    %plot(t(1:n),H(1:n));
    %plot(t(1:n),Error1(1:n));
    %plot(t(1:n),Error2(1:n));
end