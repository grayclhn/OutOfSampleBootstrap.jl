module wip
using OutOfSampleBootstraps
using Docile

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
