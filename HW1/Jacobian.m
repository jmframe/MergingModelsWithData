function [J] = Jacobian(X, M, step, m , n)
    Mnew = M * (1+step);
    J = zeros(m,n);
    for i = 1:m
    for j = 1:n
        J(i,j) = (X(i,j) * Mnew(j,1) - X(i,j) * M(j,1)) / (step*M(j,1));
        %pause(0.05)
    end
    end
end