module DCATools2

using Random

import Base: size

export PottsGraph
export rand_pottsgraph

include("objects.jl")
include("sample.jl")

end # module
