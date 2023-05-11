function [y] = RK4_with_step_control(dx1,dx2,x1,x2,h,a,b,epsW,epsB,hmin,imax)
% Returns y, which represents a list of points (x1,x2) calculated using
% ADAMS PREDICTOR-CORRECTOR method based on provided differential equations
% (dx1,dx2) for t c [a,b]
%  
%   INPUT
%       dx1     -  first differential equation
%       dx2     -  seconds differential equation
%       x1      -  initial x1 (dx1(a))
%       x2      -  initial x2 (dx2(a))
%       h       -  step lenght
%       a       -  start point
%       b       -  stop point
%       epsW    -  relative measurement 
%       epsB    -  absolute measurement error
%       hmin    -  minimal step lenght
%       imax    -  limit of iterations
%
%   OUTPUT
%       y       -  vector of points (x1,x2) 
%
%SAMPLE RUN 
%       >> [y] = RK4_with_step_control(@(x1,x2)x2 + x1*(0.5-(x1)^2-(x2)^2),@(x1,x2)-x1 + x2*(0.5-(x1)^2-(x2)^2),0.03,0.001,0.1,0,20,1e-6,1e-6,0.1,200)
%
%
    tic
    t = zeros(imax,1); % time [a,b]
    t(1,1) = a;
    X1=zeros(imax,1);
    X2=zeros(imax,1);
    Error1=zeros(imax,1);
    Error2=zeros(imax,1);
    X1(1,1)=x1;
    X2(1,1)=x2;
    x1d=x1;
    x2d=x2;
    H=zeros(imax,1); % steps length H(t)
    H(1,1)=h;
    n=1;
    while(t(n,1)<b && n < imax)
        [x1,x2]=rk4(dx1,dx2,x1d,x2d,h); % calculating using rk4 algorythm (k1,k2,k3,k4)

        [x1d,x2d]=rk4(dx1,dx2,x1d,x2d,h/2);
        [x1d,x2d]=rk4(dx1,dx2,x1d,x2d,h/2);
        X1(n+1,1)=x1d;
        X2(n+1,1)=x2d;

        errorX1=(x1d-x1)/15; % calculating aproximation error
        errorX2=(x2d-x2)/15;
        Error1(n+1,1)=errorX1;
        Error2(n+1,1)=errorX2;
        eps1 = abs(x1d) * epsW + epsB; % calculating measurement error
        eps2 = abs(x2d) * epsW + epsB;
    
        alfa1 = (eps1/abs(errorX1))^(1/5);
        alfa2 = (eps2/abs(errorX2))^(1/5);
        alfa = min(alfa1, alfa2);
        if (0.9 * alfa >= 1)
            if (t(n,1) + h < b)
                corr = 0.9 * alfa * h;
                h = min(b - t(n,1), min(5*h, corr));
            end
        else
            corr = 0.9 * alfa * h;
            if(corr >= hmin)
                h = corr;
            end
        end
        t(n+1,1) = t(n,1)+h;
        H(n+1,1)=h;
        n=n+1;
    end
    toc;
    y = [X1(1:n), X2(1:n)];
    showPlot(X1, X2, t, H, n, Error1,Error2);
end

function [x1,x2] = rk4(dx1,dx2,x1,x2,h)
    halfstep=h/2;

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
end


function [] = showPlot(X1, X2, t, H, n, Error1,Error2)
    hold on
    xlabel('x1');
    ylabel('x2');
    plot(X1(1:n),X2(1:n));
    %xlim([0, 20]);
    %plot(t(1:n),X1(1:n));
    %plot(t(1:n),X2(1:n));
    %plot(t(1:n),H(1:n));
    %plot(t(1:n),Error1(1:n));
    %plot(t(1:n),Error2(1:n));
    title('RK4 with step control')
end