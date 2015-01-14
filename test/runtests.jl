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
