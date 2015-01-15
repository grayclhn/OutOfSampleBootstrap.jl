using OutOfSampleBootstraps
using Base.Test

n = 100
R = 30
P = n - R
k = 4

y = randn(n)
x = [ones(n) randn(n, k-1)]

estimates = Array(Float64, k, P)
errors = Array(Float64, P)

## Verify that recursive_ols produces the right estimates.
b1 = recursive_ols(y, x, R)
@test b1[:,1] == x[1:R,:] \ y[1:R]
@test b1[:,end] == x[1:end-1,:] \ y[1:end-1]

## Make sure that different functions to get the recursive OLS
## estimates return the same values.
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
