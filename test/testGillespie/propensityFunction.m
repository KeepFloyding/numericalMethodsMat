function a = propensityFunction(X)

%==================================================================
% GENERATING PROPENSITY FUNCTION
%==================================================================

global k_f k_r


a(1) = k_f*X(1)*X(2);
a(2) = k_r*X(3);

end