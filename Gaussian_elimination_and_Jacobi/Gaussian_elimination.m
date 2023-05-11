function X = Task_1()
    n=[5,10,25,50,100,200];
    approaches=size(n, 2);
    x=n;
    yA=zeros(1,approaches);
    yB=zeros(1,approaches);
    T = zeros(1, approaches);
    for i=1:approaches
        [A, b]=matrix_a(n(i));
        fprintf('Time of solving %d equations:\n', n(i));
        tic
        [X, epsilon] = Gauss_partial_pivoting(A,b, n(i));
        T(i)=toc;
        fprintf('Measurement error: %d\n\n', epsilon);
        yA(i)=epsilon;
    end
    for i=1:approaches
        [A, b]=matrix_b(n(i));
        fprintf('Time of solving %d equations:\n', n(i));
        [X, epsilon] = Gauss_partial_pivoting(A,b);
        fprintf('Measurement error: %d\n\n', epsilon);
        yB(i)=epsilon;
    end
    figure(1);
    plot(x, yA);
    title('Measurement error for n equations (A)');
    xlabel('Number of data (n)');
    ylabel('Measurement error');
    
    figure(2);
    plot(x, yB)
    title('Measurement error for n equations (B)');
    xlabel('Number of data (n)');
    ylabel('Measurement error');

    figure(3);
    plot(x, T);
    title('Time of solving n equations (A)');
    xlabel('Number of data (n)');
    ylabel('Time of solving');

    clear X;
end


function [X, epsilon] = Gauss_partial_pivoting(A, b, n)
tic
    if nargin < 3
        n = size(A, 1);
    end
    for i=1:n
        maxA=abs(A(i,i)); % sets initial maximum value i column
        maxi=i;
        for j=i+1:n
            if abs(A(j,i))>maxA
                maxA=abs(A(j,i));
                maxi=j;
            end
        end
        temp = A(i,:); % switching rows
        A(i,:) = A(maxi,:);
        A(maxi,:) = temp;
        tempB = b(i); % switching rows
        b(i) = b(maxi);
        b(maxi) = tempB;
        for j=i+1:n
            l = A(j,i)/A(i,i);
            A(j,:) = A(j,:)-l*A(i,:);
            b(j) = b(j) - l*b(i);
        end
    end
    X = zeros(n,1);
    for i=n:-1:1
        equation=A(i,:)*X;
        X(i,1) = (b(i,1)-equation)/A(i,i);
    end
    z=A*X-b;
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

function [A, b] = matrix_b(n)
    A=zeros(n,n);
    b=zeros(n,1);
    for i=1:n
        for j=1:n
            if i == j
                A(i, j) = 1/6;
            else
                A(i, j) = 2*(i-j)+2;
            end
        end
        b(i,1) = 2.5 + 0.4 * i;
    end
end
