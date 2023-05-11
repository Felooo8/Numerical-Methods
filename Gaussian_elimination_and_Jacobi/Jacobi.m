function X = Task_2()
    n=[5,10,25,50,100,200];
    approaches=size(n, 2);
    x=n;
    yA=zeros(1,approaches);
    T = zeros(1, approaches);
    for i=1:approaches
        [A, b]=matrix_a(n(i));
        tic
        fprintf('Time of solving %d equations:\n', n(i));
        [X, epsilon] = Jacobi(A,b);
        T(i)=toc;
        fprintf('Measurement error: %d\n\n', epsilon);
        yA(i)=epsilon;
    end

    figure(1);
    plot(x, yA);
    title('Measurement error for n equations (A)');
    xlabel('Number of data (n)');
    ylabel('Measurement error');
    
    figure(2);
    plot(x, T);
    title('Time of solving n equations (A)');
    xlabel('Number of data (n)');
    ylabel('Time of solving');

    clear X
end


function [X, epsilon] = Jacobi(A, b, epsilon, n)
    tic
    if nargin < 4
        n = size(A, 1);
        if nargin < 3
            epsilon = 10^(-6);
        end
    end
    % we will only need inverse of D martix
    invD = zeros(n);
    % we will only need the sum of L and U matrixes
    LU = zeros(n);
    for i=1:n
        for j=1:n
            if j~=i
                LU(i,j) = A(i,j);
            end
        end
        % D is diagonal matrix, so to inverse it we only have to inverse numbers
        invD(i,i) = 1/A(i,i);
    end
    Xj=zeros(n,1);
    euklides = 1;
    
    while euklides > epsilon
        Xi=Xj;
        Xj = -invD*LU*Xi+invD*b;
        euklides = norm(Xj-Xi);
    end
    X=Xj;
    epsilon=norm(A*X-b);
    toc
end


function [A, b] = matrix_a(n)
    A=zeros(n,n);
    b=zeros(n,1);
    for i=1:n
        for j=1:n
            if i == j
                A(i, j) = 10;
            elseif j == i - 3 || j == i + 3
                A(i, j) = -3+j/n;
            end
        end
        b(i,1) = 2.5+0.5 * i;
    end
end
