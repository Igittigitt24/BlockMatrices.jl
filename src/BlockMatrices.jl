module BlockMatrices

export ZeroMatrix, IdentityMatrix, BlockMatrix

abstract type AbstractBlockMatrix{T} <: AbstractArray{T,2} end

Base.size(A::AbstractBlockMatrix) = (A.m, A.n)

struct ZeroMatrix{T} <: AbstractMatrix{T}
    m::Int
    n::Int
end
Base.size(A::ZeroMatrix) = (A.m, A.n)

function Base.A_mul_B!{T}(Y::AbstractVector{T}, A::ZeroMatrix{T}, B::AbstractVector{T})
    @assert size(A,2) == size(B,1)
    @assert size(A,1) == size(Y,1)
    fill!(Y, zero(T))
end
function Base.At_mul_B!{T}(Y::AbstractVector{T}, A::ZeroMatrix{T}, B::AbstractVector{T})
    @assert size(A,1) == size(B,1)
    @assert size(A,2) == size(Y,1)
    fill!(Y, zero(T))
end
Base.Ac_mul_B!{T}(Y::AbstractVector{T}, A::ZeroMatrix{T}, B::AbstractVector{T}) = At_mul_B!(Y,A,B)

struct IdentityMatrix{T} <: AbstractMatrix{T}
    m::Int
    λ::T
    IdentityMatrix{T}(λ::T,s) where {T} = new(s,λ)
end
Base.size(A::IdentityMatrix) = (A.m, A.m)

function Base.A_mul_B!{T}(Y::AbstractVector{T}, A::IdentityMatrix{T}, B::AbstractVector{T})
    @assert size(A,2) == size(B,1)
    @assert size(A,1) == size(Y,1)
    copy!(Y, B)
    scale!(Y, A.λ)
    return Y
end
Base.At_mul_B!{T}(Y::AbstractVector{T}, A::IdentityMatrix{T}, B::AbstractVector{T}) = A_mul_B!(Y, A, B)
Base.Ac_mul_B!{T}(Y::AbstractVector{T}, A::IdentityMatrix{T}, B::AbstractVector{T}) = A_mul_B!(Y, A, B)

"""
    M = BlockMatrix{T}(Mtuple) <: AbstractBlockMatrix{T} <: AbstractArray{T,2}

BlockMatrix of form `[A11 A12 ...; A21 A22 ...; ...]` where `Mtuple = ((A11,A12,...),(A21,A22,...),...)`
and each `Aij<:AbstractArray{T}``
"""
struct BlockMatrix{T, MT<:NTuple{<:Any,NTuple{<:Any,AbstractArray}}} <: AbstractBlockMatrix{T}
    m::Int
    n::Int
    M::MT
    mranges::Array{UnitRange{Int},1}
    nranges::Array{UnitRange{Int},1}
    tmp::Array{T,1}         #Size: max(m,n)
    function BlockMatrix{T,MT}(M::MT) where {T, mi, ni, MT <: NTuple{mi,NTuple{ni,AbstractArray{T}}}}
        hs = Array{Int,1}(mi)
        ws = Array{Int,1}(ni)
        for i = 1:mi                #All rows have equal hight
            h = size(M[i][1],1)
            for j = 1:ni
                @assert size(M[i][j],1) == h
            end
            hs[i] = h
        end
        for j = 1:ni                #All columns have equal width
            w = size(M[1][j],2)
            for i = 1:mi
                @assert size(M[i][j],2) == w
            end
            ws[j] = w
        end
        cumh = cumsum([0;hs])
        cumw = cumsum([0;ws])
        mranges = [(cumh[i]+1):cumh[i+1] for i = 1:mi]
        nranges = [(cumw[j]+1):cumw[j+1] for j = 1:ni]
        new(sum(hs), sum(ws), M, mranges, nranges, Array{T,1}(max(sum(hs),sum(ws))))
    end
end
BlockMatrix(M::MT) where {T, mi, ni, MT <: NTuple{mi,NTuple{ni,AbstractArray{T}}}} = BlockMatrix{T,MT}(M)

Base.showarray(io::IO, X::BlockMatrix, b::Bool) = display(io, X)
function Base.display(io::IO, M::BlockMatrix{T,MT}) where {T, mi, ni, MT <: NTuple{mi,NTuple{ni,AbstractArray{T}}}}
    println(io, "$(size(M,1))x$(size(M,2)) BlockMatrix{$T}:")
    println(io, "With blocks of types:")
    Base.print_matrix(io, vcat((hcat(typeof.(Mi)...) for Mi in M.M)...), "[", ", ", "]")
    println(io, "\nand sizes:")
    Base.print_matrix(io, vcat((hcat(size.(Mi)...) for Mi in M.M)...), "[", ", ", "]")
end

import Base.Cartesian: @nexprs

@generated function Base.A_mul_B!(Y::AbstractVector{T}, M::BlockMatrix{T, MT}, B::AbstractVector{T}) where {T, mi, ni, MT <: NTuple{mi,NTuple{ni,AbstractArray{T}}}}
    ni1 = ni-1
    quote
        @assert size(B,2) == size(Y,2) == 1
        @assert size(M,2) == size(B,1)
        @assert size(M,1) == size(Y,1)

        @nexprs $mi i -> A_mul_B!(view(Y, M.mranges[i]), M.M[i][1], view(B,M.nranges[1]))

        @nexprs $ni1 j -> begin
            @nexprs $mi i -> A_mul_B!(view(M.tmp, M.mranges[i]), M.M[i][j+1], view(B,M.nranges[j+1]))
             Y .+= view(M.tmp, 1:M.m)
        end
        return Y
    end
end

@generated function Base.At_mul_B!(Y::AbstractVector{T}, M::BlockMatrix{T, MT}, B::AbstractVector{T}) where {T, mi, ni, MT <: NTuple{mi,NTuple{ni,AbstractArray{T}}}}
    mi1 = mi-1
    quote
        @assert size(B,2) == size(Y,2) == 1
        @assert size(M,1) == size(B,1)
        @assert size(M,2) == size(Y,1)

        @nexprs $ni j -> At_mul_B!(view(Y, M.nranges[j]), M.M[1][j], view(B,M.mranges[1]))

        @nexprs $mi1 i -> begin
            @nexprs $ni j -> At_mul_B!(view(M.tmp, M.nranges[j]), M.M[i+1][j], view(B,M.mranges[i+1]))
             Y .+= view(M.tmp, 1:M.n)
        end
        return Y
    end
end


@generated function Base.Ac_mul_B!(Y::AbstractVector{T}, M::BlockMatrix{T, MT}, B::AbstractVector{T}) where {T, mi, ni, MT <: NTuple{mi,NTuple{ni,AbstractArray{T}}}}
    mi1 = mi-1
    quote
        @assert size(B,2) == size(Y,2) == 1
        @assert size(M,1) == size(B,1)
        @assert size(M,2) == size(Y,1)

        @nexprs $ni j -> Ac_mul_B!(view(Y, M.nranges[j]), M.M[1][j], view(B,M.mranges[1]))

        @nexprs $mi1 i -> begin
            @nexprs $ni j -> Ac_mul_B!(view(M.tmp, M.nranges[j]), M.M[i+1][j], view(B,M.mranges[i+1]))
             Y .+= view(M.tmp, 1:M.n)
        end
        return Y
    end
end

Base.size(A::BlockMatrix) = (A.m, A.n)

end # module
