using IsoTNS, Test
using IsoTNS: @ein_str
using LinearAlgebra

@testset "First test set" begin
    phix = 2
    a = 2
    b = 2
    c = 2
    i = c

    nSTT1 = rand(ComplexF64, phix, a, 1, b, c)
    nSTT2 = rand(ComplexF64, phix, a, i, b, c)
    nSTT3 = rand(ComplexF64, phix, a, i, b, 1)
    nSTT = [nSTT1, nSTT2, nSTT3]


    MM1, MM2 = MM_line!(MMPS(nSTT), 100, 0)
    @test size(MM1.tensors[1]) == (a, a, 1, b, c)

    pp = ein"rqnsi, rqnti ->st "(MM1.tensors[3], conj.(MM1.tensors[3]))
    
    identity_matrix = Matrix{Float64}(I, 4, 4)
    @test pp ≈ identity_matrix

end