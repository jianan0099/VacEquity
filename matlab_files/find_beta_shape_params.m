function [a,b] = find_beta_shape_params(mu,sd)
  var = sd ^ 2;
  a = ((1 - mu) / var - 1 / mu) * mu ^ 2;
  b = a * (1 / mu - 1);
end

