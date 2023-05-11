function [y] = P_K_Adams(dx1,dx2,x1,x2,h,a,b)
% Returns y, which represents a list of points (x1,x2) calculated using
% ADAMS PREDICTOR-CORRECTOR method based on provided differential equations
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
%       >> [y] = P_K_Adams(@(x1,x2)x2 + x1*(0.5-(x1)^2-(x2)^2),@(x1,x2)-x1 + x2*(0.5-(x1)^2-(x2)^2),0.03,0.001,0.01,0,20)
%
%
    tic;
    t=a:h:b;
    X1(:,1) = x1;
    X2(:,1) = x2;
    % initializing
    for i=1:3
        [x1,x2]=rk4(dx1,dx2,x1,x2,h); % calculating using rk4 algorythm (k1,k2,k3,k4)

        X1(:,i+1)=x1;
        X2(:,i+1)=x2;
    end
    % constants
    B_bashforth = [55/24,-59/24,37/24,-9/24];
    B_moulton = [252/720,646/720,-264/720,106/720,-19/720];
    for i = 4:(length(t)-1)
        % prediction
        sum1 = 0;
        sum2 = 0;
        for j=0:3
            sum1 = sum1 + B_bashforth(j+1)*dx1(X1(i-j),X2(i-j));
            sum2 = sum2 + B_bashforth(j+1)*dx2(X1(i-j),X2(i-j));
        end
        pred1 = x1 + sum1*h;
        pred2 = x2 + sum2*h;
        % corection
        sum1 = 0;
        sum2 = 0;
        for j=0:3
            sum1 = sum1 + B_moulton(j+2)*dx1(X1(i-j),X2(i-j));
            sum2 = sum2 + B_moulton(j+2)*dx2(X1(i-j),X2(i-j));
        end
        x1 = x1 + (sum1+B_moulton(1)*dx1(pred1, pred2))*h;
        x2 = x2 + (sum2+B_moulton(1)*dx2(pred1, pred2))*h;

        X1(:,i+1)=x1;
        X2(:,i+1)=x2;
    end
    toc;
    y = [X1; X2];
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