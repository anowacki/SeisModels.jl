using Test

@testset "All tests" begin
    include("models_1d.jl")
    include("inbuilt_models.jl")
    include("basic_properties.jl")
end
