module SpatialRBFEstimation

using LinearAlgebra, GeometryBasics, Distances, RecipesBase

struct RadialBasisFunction

    weights::Vector{T} where T<:Number
    nodes::Vector{Point2{T}} where T<:Number
    bandwidths::Vector{T} where T<:Number


    function RadialBasisFunction(weights::Vector{T}, nodes::Vector{Point2{T}}, bandwidths::Vector{T}) where T<:Number
        @assert length(weights) == length(nodes) == length(bandwidths)
        @assert minimum(weights)>0 "Weights must be positive"
        @assert minimum(bandwidths)>0 "Bandwidths must be positive"
        new(weights,nodes,bandwidths)
    end

    function RadialBasisFunction(weight::T, node::Point2{T}, bandwidth::T) where T<:Number
        RadialBasisFunction([weight],[node],[bandwidth])
    end

    function RadialBasisFunction()
        RadialBasisFunction(1.0,Point2{Float64}(0.0,0.0),1.0)
    end

end

function Base.:+(functions::RadialBasisFunction...)

    weigths = reduce(vcat,[f.weights for f in functions])
    nodes = reduce(vcat,[f.nodes for f in functions])
    bandwidths = reduce(vcat,[f.bandwidths for f in functions])

    return RadialBasisFunction(weigths,nodes,bandwidths)
end

function Base.:*(k::Number, f::RadialBasisFunction)

    return RadialBasisFunction(k*f.weights,f.nodes,f.bandwidths)
end

function (f::RadialBasisFunction)(x::Point)
    
    w = f.weights
    ker = exp.(-Euclidean().(f.nodes,Ref(x)).^2 ./(2*f.bandwidths))
    return sum(w.*ker)

end

function (f::RadialBasisFunction)(x::Vector{T}) where T<:Number
    
    p = Point2(x)
    return f(p)

end

@recipe function f(::Type{RadialBasisFunction},rbf::RadialBasisFunction)

    x = collect(-1:.01:1)
    y = collect(-1:.01:1)
    val = [rbf([xx,yy]) for xx in x, yy in y]

    title --> "RBF density"
    aspectratio --> :equal
    colorbar --> false
    @series begin
        x,y,val
    end

end


end
