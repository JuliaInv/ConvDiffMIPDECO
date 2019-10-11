export dilation

"""
function dilation(mc,M::AbstractMesh,ssize=2)

dilated the binary-valued image mc by ssize cells.
The code works in 2D and 3D

"""
function dilation(mc::BitArray,M::AbstractMesh,ssize=2,boundLow=0.1,boundHigh=0.9)

    n   = M.n
    mc  = reshape(copy(mc),tuple(n...))

    A   = ones(1+2*ssize)
    A /= sum(vec(A))
    rId = Array{Float64}(mc)

    perm = zeros(Int64,M.dim)
    perm[1:M.dim-1] = collect(2:M.dim)
    perm[end]=1
    perm
    for k=1:M.dim
        szrId = size(rId)
        rId = reshape(rId,size(rId,1),:)
        for j=1:size(rId,2)
            rId[:,j] = conv(A,rId[:,j])[1+ssize:end-ssize]
        end
        rId = reshape(rId,szrId)
        rId = permutedims(rId,perm)
    end
    rId = (vec(rId).>=boundLow) .& (vec(rId) .<= boundHigh)
    return rId
end
