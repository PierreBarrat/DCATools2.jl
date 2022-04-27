struct PottsGraph
	J :: Array{Float64, 4} 
	h :: Array{Float64, 2} 
end

PottsGraph(L::Int, q::Int) = PottsGraph(zeros(Float64, q, q, L, L), zeros(Float64, q, L))

size(g::PottsGraph) = size(g.h)