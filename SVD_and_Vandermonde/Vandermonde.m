function y = task_2(x, Xdata, Ydata, showPlot)
% Returns y=f(x) for given x and optinally shows plot(if showPlot is true, default=false)
% using Vandermonde decomposition and interpolation
% Example of parameters:
% x=-4, 
% Xdata= [-10,-8,-6,-4,-2,0,2,4,6,8,10], 
% Ydata=[-18.7370;-8.1583;-1.9146;-0.3887;1.8030;1.1890;0.4738; 0.4726;
%           0.0941;-2.3716;-6.6512],
% showPlot =  true

    if nargin < 4
        showPlot = false;
    end

    X = fliplr(vander(Xdata)); % create Vandermonde's matrix + flip
    a=inv(X)*Ydata;
    s=length(Xdata);
    
    % calculating y=f(x)
    y = a(1); % y = absolute term
    for i=2:s
        y = y + a(i)*x.^(i-1); % adding each monomial
    end


    if showPlot % shows the plot of polynomial + measurement error
        show_plot(Xdata, Ydata, a, X, s);
    end
    
end


function [] = show_plot(Xdata, Ydata, a, X, s)
    samp = 20;
    sampling = s*samp;
    x=linspace(min(Xdata)-0.5,max(Xdata)+0.5, sampling);
    Y=task_2(x, Xdata, Ydata, false); % matrix of y=f(x) for each x    
    
    hold on
    plot(Xdata, Ydata, 'ob', 'DisplayName', 'Data');
    title('Polynomial');
    xlim([-11, 11]);
    plot(x, Y, 'g-', 'LineWidth', 2, 'DisplayName', 'Polynomial');
    legend('Location', 'southeast');
    
    [euk, maxi] = return_epsilon(a, X, Ydata);

    fprintf('Euklides measurement error for polynomial: %d\n', euk);
    fprintf('Maximum measurement error for polynomial: %d\n', maxi);

end

function [euk, maxi] = return_epsilon(A, X, Ydata)
    res = Ydata-X*A;
    euk = norm(res, 2);
    maxi = norm(res, Inf);
end