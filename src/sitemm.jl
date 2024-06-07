struct NMPS{T}
    tensors::Vector{Array{T, 6}}
end

function truncated_svd(tmat::AbstractMatrix, dmax::Int, atol::Real)
    u, s, v = LinearAlgebra.svd(tmat)
    dmax = min(searchsortedfirst(s, atol, rev=true), dmax, length(s))
    return u[:, 1:dmax], s[1:dmax], v'[1:dmax, :], sum(s[dmax+1:end])
end

function sperate_singleT_tripartite(site::NMPS,dmax::Int, atol::Real)
    nSTT = site.tensors[1]
    re_nSTT = reshape(nSTT, size(nSTT, 1)*size(nSTT, 2)*size(nSTT, 3), :)
    u1, s1, v1 = truncated_svd(re_nSTT, dmax, atol)
    A = reshape(u1, size(nSTT, 1), size(nSTT, 2), size(nSTT, 3), size(nSTT, 2), :)

    nSTT_B_2 = reshape(Diagonal(s1) * v1,  size(nSTT, 2), :, size(nSTT, 4), size(nSTT, 5), size(nSTT, 6))
    nSTT_B_2 = permutedims(nSTT_B_2, (2, 3, 4, 5, 1))
    re_nSTT_B_2 = reshape(nSTT_B_2, size(nSTT_B_2, 1)*size(nSTT_B_2, 2)*size(nSTT_B_2, 3), :)
    u2, s2, v2 = truncated_svd(re_nSTT_B_2, dmax, atol)

    B = reshape(u2,  :, size(nSTT, 4), size(nSTT,5), size(s2, 1))
    C = reshape(Diagonal(s2) * v2, size(s2, 1), size(nSTT, 6), :)

    return A, B, C
end