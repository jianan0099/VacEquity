function sample_results = get_trunc_beta(mu,sd,low,high)

[a,b] = find_beta_shape_params(mu, sd);
pd= truncate(makedist('Beta',a,b),low,high);
sample_results = random(pd,100000,1);
end

