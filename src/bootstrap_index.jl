export Block

type Block{T}
    length::Integer
end

Block(s::Symbol, n) = Block{s}(n)

@doc """bootstrap_index! produces an array of integers that produces a draw from
the bootstrap distribution when used to index the original data. The first
argument is a vector of length n or an n by m matrix that will store the index.
The integer n is the number of observations and m is the number of bootstrap
replications. The second argument is a `Block` object that specifies the
blocking method and the block length.

Only the circular and moving blocks bootstrap are supported, which can be
invoked by using `Block(:circular, n)` and `Block(:moving, n)` respectively.

Examples
--------
i = Array(Int, 90) ## initialization

bootstrap_index!(i, Block(:circular, 10))
bootstrap_index!(i, Block(:moving, 9)) ## Overwrites previous line's results
""" ->
function bootstrap_index!{T<:Integer}(index::AbstractVector{T}, block::Block{:circular})
    n = length(index)
    nblocks = cld(n, block.length)
    for i in 1:nblocks
        blockstart_i = rand(1:n)
        for t in ((i-1) * block.length + 1):min(i*block.length, n)
            ## We're being sort of sloppy with the offsets here. The
            ## key issue is that there is a single startpoint for each
            ## block, and that starting point is uniformly distributed
            ## on 1:n. Even though t increases as i increases,
            ## randomizing blockstart_i ensures this uniform
            ## distribution.
            index[t] = mod(blockstart_i + t, n) + 1
        end
    end
end

function bootstrap_index!{T<:Integer}(index::AbstractVector{T}, block::Block{:moving})
    n = length(index)
    nboot = 0
    while nboot < n
        blockstart_i = rand(0:n-1)
        ## Shorten the block if it is larger than the remaining spots
        ## in index or if it runs over the end of the dataset.
        blocklength_i = min(n - blockstart_i, n - nboot, block.length)
        index[nboot + (1:blocklength_i)] = blockstart_i + (1:blocklength_i)
        nboot += blocklength_i
    end
end

function bootstrap_index!{T<:Integer}(index::AbstractMatrix{T}, block::Block)
    for i in 1:size(index, 2)
        bootstrap_index!(sub(index,:,i), block)
    end
end
