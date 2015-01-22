export recursive_ols, recursive_ols!

recursive_ols_overview =
"""Calculates the OLS estimators for the model y ~ x recursively and returns
them as an k by n - R matrix. The jth column of `estimates` is the OLS
estimator estimated by regressing y[1:(R+j-1)] on x[1:(R+j-1),:] and is used to
forecast y[R+j] rom x[R+j,:].

* `recursive_ols` creates a new matrix containing the estimates and does basic
  checks that the arrays have the right dimensions. The actual calculations are
  done by `recursive_ols!`

* `recursive_ols!` fills its first argument with the OLS estimates and does not
  check the dimensions of its arguments. It is intended to be used inside
  simulations, etc., where the arguments are known to be safe.

There are also a methods that return the OLS estimators and the implied
prediction errors.
"""
Index = Union(AbstractVector{Int},Range{Int},UnitRange{Int})

@doc recursive_ols_overview ->
function recursive_ols!(estimates::AbstractMatrix{Float64},
                        y::AbstractVector{Float64},
                        x::AbstractMatrix{Float64},
                        order::Index)
    k, P = size(estimates)
    R = length(y) - P
    ## QR decomposition used to construct the (hypothetical) period 0
    ## estimate
    Rmat = vcat(qrfact(hcat(x[order[1:R-1],:], y[order[1:R-1]]))[:R], zeros(1, k+1))
    for j in 1:P
        ols_update!(sub(estimates,:,j), Rmat, x, y, order[j+R-1])
    end
end

function recursive_ols!(estimates::AbstractMatrix{Float64},
                        errors::AbstractVector{Float64},
                        y::AbstractVector,
                        x::AbstractMatrix,
                        order::Index)
    k, P = size(estimates)
    R = length(y) - P
    ## QR decomposition used to construct the (hypothetical) period 0
    ## estimate
    Rmat = vcat(qrfact(hcat(x[order[1:R-1],:], y[order[1:R-1]]))[:R], zeros(1, k+1))
    for j in 1:P
        ols_update!(sub(estimates,:,j), Rmat, x, y, order[j+R-1])
        errors[j] = y[order[R+j]] - scalar(x[order[R+j],:] * estimates[:,j])
    end
end    

@doc recursive_ols_overview ->
function recursive_ols(y::AbstractVector, x::AbstractMatrix, R::Integer, order::Index)
    length(y) == size(x,1) ||
      error("`y` should be a vector with the same number of elements as rows of `x`")
    estimates = Array(Float64, size(x, 2), length(y) - R)
    recursive_ols!(estimates, y, x, order)
    return estimates
end 

@doc "Function that actually updates the OLS coefficients." ->
function ols_update!(b::AbstractVector{Float64}, R::AbstractMatrix{Float64},
                     x::AbstractMatrix{Float64}, y::AbstractVector{Float64},
                     j::Integer)
    k = size(x, 2)
    R[k+2,1:k] = x[j,:]
    R[k+2,k+1] = y[j]
    for i in 1:k
        Gi, _ = givens(R, i, k+2, i)
        A_mul_B!(Gi, R)
    end
    b[1:k] = copy(R[1:k,k+1])
    Base.LinAlg.LAPACK.trtrs!('U', 'N', 'N', sub(R, 1:k, 1:k), b)
end

## Additional convenience methods
recursive_ols!(estimates, y, x) = recursive_ols!(estimates, y, x, 1:length(y))
recursive_ols!(estimates, errors, y, x) =
    recursive_ols!(estimates, errors, y, x, 1:length(y))
recursive_ols(y, x, R) = recursive_ols(y, x, R, 1:length(y))
