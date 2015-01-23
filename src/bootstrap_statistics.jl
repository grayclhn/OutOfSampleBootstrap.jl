export bootstrap_location, bootstrap_interval

function bootstrap_location(::OutOfSampleBootstrap{:naive},
                            ::Block,
                            f::AbstractVector)
    mean(f)
end

function bootstrap_location(::OutOfSampleBootstrap{:cs07_ols},
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
    return fsum / (P * (n-1))
end

function bootstrap_location(::OutOfSampleBootstrap{:mine_ols},
                            ::Block,
                            y::AbstractVector,
                            x::AbstractMatrix,
                            L::Function)
    βhat = x \ y
    fsum = 0.0
    for t in 1:length(y)
        fsum += L(y[t] - scalar(x[t,:] * βhat))
    end
    fsum / length(y)
end

function bootstrap_interval(α::Float64,
                            nboot::Integer,
                            boot::OutOfSampleBootstrap{:naive},
                            block::Block,
                            f::AbstractVector)
    b = Array(Float64, nboot)
    recursive_bootstrap!(b, boot, block, f)
    quantile(b, [α/2, 1 - α/2]) - bootstrap_location(boot, block, f)
end

function bootstrap_interval(α::Float64,
                            nboot::Integer,
                            boot::OutOfSampleBootstrap{:cs07_ols},
                            block::Block,
                            y::AbstractVector,
                            x::AbstractMatrix,
                            βhat::AbstractMatrix,
                            L::Function)
    b = Array(Float64, nboot)
    recursive_bootstrap!(b, boot, block, y, x, βhat, L)
    return quantile(b, [α/2, 1 - α/2]) - bootstrap_location(boot, y, x, βhat, L)
end

function bootstrap_interval(α::Float64,
                            nboot::Integer,
                            boot::OutOfSampleBootstrap{:mine_ols},
                            block::Block,
                            y::AbstractVector,
                            x::AbstractMatrix,
                            R::Integer,
                            L::Function)
    b = Array(Float64, nboot)
    recursive_bootstrap!(b, boot, block, y, x, R, L)
    return quantile(b, [α/2, 1 - α/2]) - bootstrap_location(boot, block, y, x, L)
end
