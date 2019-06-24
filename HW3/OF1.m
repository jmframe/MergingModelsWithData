function [cost] = OF1(theta)
    cost = theta(1,1) + (theta(1,2) - 1)^2;
end
