function result = g_lambda_func( z, rho, this_lambda )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
this_gamma = sqrt((2 / rho) + (1 / this_lambda));
this_kappa = 1 / this_lambda / this_gamma;
if (z > this_gamma)
    result = 1;
elseif (z <= this_lambda && z >= this_kappa)
    result = 1 - 0.5 * rho * (z - this_gamma) * (z - this_gamma);
elseif (z <= this_kappa && z >= 0)
    result = this_lambda * z * z;
else
    result = 0;
end
end

