n = 15_000
R =  1_000
k =    400
P = n - R

y = randn(n)
x = [ones(n) randn(n, k-1)]

@time b1 = recursive_ols(y, x, R)

