using IsoTNS, Test

# Ref: https://github.com/CodingThrust/SimpleTDVP.jl/blob/main/test/mps.jl
@testset "First test set" begin
    phix = 7

    a = 4
    i1 = 2
    i2 = 3
    b = 5
    c = 6

    nSTT = rand(ComplexF64, phix, a, i1, i2, b, c)

    A, B, C = sperate_singleT_tripartite(NMPS([nSTT]),min(phix, a, b, c), 0)

    @show pp
    @show nSTT

    (pp-nSTT)./nSTT
end