function y = task_1(x, degree, Xdata, Ydata, showPlot)
% Returns y=f(x) for given x, degree of polynomial, Xdata, Ydata and
% optionaly shows plot and measurement error (if showPlot is true, default=false)
% Example of parameters:
% x=-4, 
% degree=7, 
% Xdata= [-10,-8,-6,-4,-2,0,2,4,6,8,10], 
% Ydata=[-18.7370;-8.1583;-1.9146;-0.3887;1.8030;1.1890;0.4738; 0.4726;
%           0.0941;-2.3716;-6.6512],
% showPlot =  true

    if nargin < 5
        showPlot = false;
    end

    X = create_x_matrix(degree, Xdata); % creates X matrix for a given polynomial degree
    z = lsq_svd(X, Ydata); % returns A matrix
    
    % calculating y=f(x) for a given x
    a = z(end:-1:1);
    y = a(1); % y = absolute term
    for i=2:degree+1
        y = y + a(i)*x.^(i-1); % adding each monomial
    end

    if showPlot % shows the plot of polynomial + measurement error
        show_plot(Xdata, Ydata, z, X, degree);
    end
end

function Z = lsq_svd(A, yi)
    % Compute the SVD of A
    [U, S, V] = svd(A);
    % find number of nonzero singular value = rank(A)
    k = length(find(S));

    Sr = S(1:k,1:k);
    Sp = inv(Sr);
    column_zeros = zeros(k, 1);
    for i =k : size(U)-1
        Sp = [Sp column_zeros];
    end
    Z = V*(Sp*U'*yi);
end

function X = create_x_matrix(degree, xi)
% creates X matrix, which represents values of x in each monomial in
% polynomial of degree = degree
% for example:
% X = create_x_matrix(3, [-2,0,2])
% X = [-8 4 2 1
%       0 0 0 1
%       8 4 2 1]
    s = size(xi, 2);
    X=zeros(s,degree+1);
    for i=1:s
        for j=1:degree+1
            X(i,j)= xi(i)^(degree+1-j);
        end
    end
end

function [] = show_plot(Xdata, Ydata, A, X, degree)
    s=length(Xdata);
    samp = 50;
    sampling = s*samp;

    x = linspace(Xdata(1)-0.5,Xdata(end)+0.5,sampling); % +/-0.5 to display whats beyond Xdata
    y = task_1(x, degree, Xdata, Ydata);

    hold on
    plot(Xdata, Ydata, 'ob', 'DisplayName', 'Data');
    xlim([-11, 11]);
    plot(x, y, 'g-', 'LineWidth', 2);

    [euk, maxi] = return_epsilon(A, X, Ydata);

    fprintf('Euklides measurement error for polynomial of degree %d: %f\n', degree, euk);
    fprintf('Maximum measurement error for polynomial of degree %d: %f\n', degree, maxi);
end

function [euk, maxi] = return_epsilon(A, X, Ydata)
    res = Ydata-X*A;
    euk = norm(res, 2);
    maxi = norm(res, Inf);
end