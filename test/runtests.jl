using OutOfSampleBootstraps
using Base.Test

# write your own tests here
n = 100
R = 30
P = n - R
k = 4

y = randn(n)
x = [ones(n) randn(n, k-1)]

estimates = Array(Float64, k, P)
errors = Array(Float64, P)

## Make sure that different functions to get the recursive OLS
## estimates return the same values.

b1 = recursive_ols(y, x, R)
recursive_ols!(estimates, y, x)
@test estimates == b1

recursive_ols!(estimates, errors, y, x)
@test estimates == b1

## Make sure basic bootstraps run
b = Array(Float64, 199) ## Initialize storage. Each call will
                        ## overwrite the previous results
recursive_bootstrap!(b, :naive, Block(:moving, 15), y, x, R)
recursive_bootstrap!(b, :cs07_ols, Block(:moving, 15), y, x, R)
recursive_bootstrap!(b, :mine_ols, Block(:moving, 15), y, x, R)
