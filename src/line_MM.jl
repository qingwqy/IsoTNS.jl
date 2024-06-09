using OMEinsum, LinearAlgebra

struct NMPS{T}
    tensors::Vector{Array{T, 6}}
end
struct MMPS{T}
    tensors::Vector{Array{T, 5}}
end

struct FMPS{T}
    tensors::Vector{Array{T, 4}}
end

struct OMPS{T}
    tensors::Vector{Array{T, 3}}
end

function truncated_svd(tmat::AbstractMatrix, dmax::Int, atol::Real)
    u, s, v = LinearAlgebra.svd(tmat)
    dmax = min(searchsortedfirst(s, atol, rev=true), dmax, length(s))
    return u[:, 1:dmax], s[1:dmax], v'[1:dmax, :], sum(s[dmax+1:end])
end

function sperate_singleT_tripartite(site::NMPS, dmax::Int, atol::Real)
    nSTT = site.tensors[1]
    re_nSTT = reshape(nSTT, size(nSTT, 1)*size(nSTT, 2)*size(nSTT, 3), :)
    u1, s1, v1 = truncated_svd(re_nSTT, dmax, atol)
    A = reshape(u1, size(nSTT, 1), size(nSTT, 2), size(nSTT, 3), :, size(nSTT, 6))

    nSTT_B_2 = reshape(Diagonal(s1) * v1,  size(nSTT, 6), :, size(nSTT, 4), size(nSTT, 5), size(nSTT, 6))
    nSTT_B_2 = permutedims(nSTT_B_2, (2, 3, 4, 5, 1))
    re_nSTT_B_2 = reshape(nSTT_B_2, size(nSTT_B_2, 1)*size(nSTT_B_2, 2)*size(nSTT_B_2, 3), :)
    u2, s2, v2 = truncated_svd(re_nSTT_B_2, dmax, atol)

    B = reshape(u2,  :, size(nSTT, 4), size(nSTT,5), size(s2, 1))
    C = reshape(Diagonal(s2) * v2, size(s2, 1), size(nSTT, 6), :)

    return A, B, C
end



function MM_line!(mpsline_O::MMPS, dmax::Int, atol::Real)
    l = length(mpsline_O.tensors)
    mpsline_G = FMPS([rand(ComplexF64, 1, 1, 1, 1) for i in 1:l])

    Imps = OMPS([rand(ComplexF64, 1, 1, 1)])
    Anx = ein"ijsmn, ksl->ijklmn"(mpsline_O.tensors[1], Imps.tensors[1])
    Anx = NMPS([Anx])
    for i in 1:l-1
        A = Anx.tensors[1]
        #A  = reshape(A, size(A, 1), size(A, 2), size(A, 3), 1, size(A, 4), size(A, 5))
        A_MM, B_MM, C_MM = sperate_singleT_tripartite(NMPS([A]), dmax, atol)
        B = mpsline_O.tensors[i+1]
        B = ein"ijsmn, ksl->ijklmn"(B,C_MM)
        mpsline_O.tensors[i] = A_MM
        Anx.tensors[1] = B
        mpsline_G.tensors[i] = B_MM
    end
    L = Anx.tensors[1]
    A_MM, B_MM, C_MM = sperate_singleT_tripartite(NMPS([L]), dmax, atol)
    mpsline_O.tensors[l] = A_MM
    mpsline_G.tensors[l] = B_MM
    return mpsline_O, mpsline_G
end

