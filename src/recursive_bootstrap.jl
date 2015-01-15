export OutOfSampleBootstrap, recursive_bootstrap!

## Simple dummy object type that will allow us to dispatch on the name
## of the bootstrap.
type OutOfSampleBootstrap{T}
end
OutOfSampleBootstrap(s::Symbol) = OutOfSampleBootstrap{s}()

## Main implementations of the different oos bootstraps. Different
## methods are invoked by dispatching on the OutOfSampleBootstrap
## abstract type (second argument).

## Naive approach used for bootstrap reality check:
function recursive_bootstrap!(b::AbstractVector,
                              ::OutOfSampleBootstrap{:naive},
                              block::Block,
                              f::AbstractVector)
    P = length(f)
    bootindex = Array(Int, P)
    fsum = sum(f)
    for i in 1:length(b)
        fsum_i = 0.0
        bootstrap_index!(bootindex, block)
        for t_i in bootindex
            fsum_i += f[t_i]
        end
        b[i] = (fsum_i - fsum) / sqrt(P)
    end
end

## Corradi-Swanson (2007) 
function recursive_bootstrap!(b::AbstractVector,
                              boot::OutOfSampleBootstrap{:cs07_ols},
                              block::Block,
                              y::AbstractVector,
                              x::AbstractMatrix,
                              βhat::AbstractMatrix,
                              L::Function)
    P = size(βhat, 2)
    n = length(y)
    R = n - P
    bootindex = Array(Int, n)
    correction = cs07_correction(boot, y, x, βhat, L)
    for i in 1:length(b)
        fsum_boot = 0.0
        bootstrap_index!(bootindex, block)
        for j in 1:P
            mAdj_j = 2 * mean((y[1:(n-1)] - x[1:(n-1),:] * βhat[:,j]) .* x[1:(n-1),:], 1)'
            # Ugly and slow; replace with the right BLAS calls laters
            βboot_ij = ((x[bootindex[1:(R+j-1)],:]'*x[bootindex[1:(R+j-1)],:]) \
                        (x[bootindex[1:(R+j-1)],:] '* y[bootindex[1:(R+j-1)]] - mAdj_j))
            fsum_boot += L(y[bootindex[R+j]] - scalar(x[bootindex[R+j],:] * βboot_ij))
        end
        b[i] = (fsum_boot - correction) / sqrt(P)
    end
end

function cs07_correction(::OutOfSampleBootstrap{:cs07_ols},
                         y::AbstractVector,
                         x::AbstractMatrix,
                         βhat::AbstractMatrix,
                         L::Function)
    P = size(βhat, 2)
    n = length(y)
    fsum = 0.0
    for j in 1:P
        for s in 2:n
            fsum += L(y[s] - scalar(x[s,:] * βhat[:,j]))
        end
    end
    return fsum / (n-1)
end

## The new one proposed by my paper
function recursive_bootstrap!(b::AbstractVector,
                              boot::OutOfSampleBootstrap{:mine_ols},
                              block::Block,
                              y::AbstractVector,
                              x::AbstractMatrix,
                              R::Integer,
                              L::Function)
    n = length(y)
    P = n - R
    bootindex = Array(Int, n)
    correction = sum(L(y - x * (x \ y)))
    βhat_i = Array(Float64, size(x, 2), P)
    errors_i = Array(Float64, P)
    
    for i in 1:length(b)
        bootstrap_index!(bootindex, block)
        ## This is a terrible approach; I should pass the index as a
        ## separate argument and reorder in recursive_ols!
        recursive_ols!(βhat_i, errors_i, y[bootindex], x[bootindex,:])
        b[i] = (sum(L(errors_i)) - correction) / sqrt(P)
    end
end

## Several variations that wrap the main versions for
## convenience.

## Allow statistic type to be specified as a symbol
recursive_bootstrap!(b, boot::Symbol, block, y, x, R, L) =
    recursive_bootstrap!(b, OutOfSampleBootstrap(boot), block, y, x, R, L)

## Default to MSE if loss function not specified
function recursive_bootstrap!(b, boot, block, y, x, R)
    L(x) = x .* x
    recursive_bootstrap!(b, boot, block, y, x, R, L)
end

## Alternative invokations for the naive oos bootstrap
function recursive_bootstrap!(b, boot::OutOfSampleBootstrap{:naive}, block, y, x, R::Integer, L::Function)
    P = length(y) - R
    βhat = Array(Float64, size(x, 2), P)
    errors = Array(Float64, P)
    recursive_ols!(βhat, errors, y, x)
    recursive_bootstrap!(b, boot, block, L(errors))
end

## Alternative invokation for the Corradi-Swanson oos bootstrap
function recursive_bootstrap!(b, boot::OutOfSampleBootstrap{:cs07_ols}, block, y, x, R::Integer, L::Function)
    βhat = recursive_ols(y, x, R)
    recursive_bootstrap!(b, boot, block, y, x, βhat, L)
end
