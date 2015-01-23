module OutOfSampleBootstraps

using Docile
@docstrings

include("utils.jl")
include("recursive_ols.jl")
include("bootstrap_index.jl")
include("recursive_bootstrap.jl")
include("bootstrap_statistics.jl")

end # module
