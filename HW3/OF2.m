function [cost] = OF2(theta)
    cost = theta(1,2) + (theta(1,1) - 1)^2;
end
