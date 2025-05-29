% Function to calculate blending functions for Bezier curves
function B = bezier_blend(t, n)
    % Calculate binomial coefficients
    C = factorial(n) ./ (factorial(0:n) .* factorial(n:-1:0));
    
    % Calculate blending functions using binomial coefficients and parameter values
    B = C .* t.^(0:n) .* (1 - t).^(n:-1:0);
end