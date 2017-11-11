function a = propensityFunction(X)

%==================================================================
% GENERATING PROPENSITY FUNCTION
%==================================================================

global k_f k_r

a = zeros(1, length(X));

a(1) = -k_f*X(1)*X(2) + k_r*X(3);
a(2) = -k_f*X(1)*X(2) + k_r*X(3);
a(3) =  k_f*X(1)*X(2) - k_r*X(3);

end