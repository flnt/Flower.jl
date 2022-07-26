struct RangeIterator
    k::Int
    d::Int
    r::Int
end

RangeIterator(n::Int, k::Int) = RangeIterator(min(n,k),divrem(n,k)...)
Base.length(it::RangeIterator) = it.k
endpos(it::RangeIterator, i::Int) = i*it.d+min(i,it.r)
Base.iterate(it::RangeIterator, i::Int=1) = i>it.k ? nothing : (endpos(it,i-1)+1:endpos(it,i), i+1)

"""
Matrix Vector multiplication
"""
function tmul2!(y::AbstractVector, A::SparseMatrixCSR, x::AbstractVector, alpha::Number, beta::Number)
    A.n == size(x, 1) || throw(DimensionMismatch())
    A.m == size(y, 1) || throw(DimensionMismatch())

    o = getoffset(A)

    @threads for row in UnitRange(1,size(y, 1))
        @inbounds begin
            accu = zero(eltype(y))
            for nz in nzrange(A, row)
                col = A.colval[nz] + o
                accu += A.nzval[nz]*x[col]
            end
            y[row] = alpha*accu + beta*y[row]
        end
    end

    return y

end

function tmul2!(y::AbstractVector, A::SparseMatrixCSR, x::AbstractVector)
    
    tmul2!(y, A, x, true, false)

end

function tmul2(A::SparseMatrixCSR, x::AbstractVector)

    y = similar(x, promote_type(eltype(A), eltype(x)), size(A, 1))

    tmul2!(y, A, x, true, false)

end

"""
Matrix Scalar multiplication
"""
function tmul2!(y::SparseMatrixCSR, A::SparseMatrixCSR, x::T, alpha::Number) where {T}
    A.n == y.n || throw(DimensionMismatch())
    A.m == y.m || throw(DimensionMismatch())

    o = getoffset(A)

    @threads for row in UnitRange(1,A.m)
        @inbounds begin
            for nz in nzrange(A, row)
                col = A.colval[nz] + o
                val = A.nzval[nz]*x
                y[row, col] = alpha*val
            end
        end
    end

    return y

end

function tmul2!(y::SparseMatrixCSR, A::SparseMatrixCSR, x::T) where {T}
    
    tmul2!(y, A, x, true)

end

function tmul2(A::SparseMatrixCSR, x::T) where {T}

    y = similar(A, promote_type(eltype(A), eltype(x)), size(A))

    tmul2!(y, A, x, true)

end

function multithread_matmul2(T::BaseThreads)

    # Matrix Vector multiplication
    @eval function  mul!(y::AbstractVector, A::SparseMatrixCSR, x::AbstractVector, alpha::Number, beta::Number)
        return tmul2!(y, A, x, alpha, beta)
    end

    @eval function  mul!(y::AbstractVector, A::SparseMatrixCSR, x::AbstractVector)
        return tmul2!(y, A, x)
    end

    @eval function  *(A::SparseMatrixCSR, x::AbstractVector)
        return tmul2(A, x)
    end

    # Matrix Scalar multiplication
    @eval function  mul!(y::AbstractVector, x::T, A::SparseMatrixCSR, alpha::Number) where {T}
        return tmul2!(y, A, x, alpha)
    end

    @eval function  mul!(y::AbstractVector, x::T, A::SparseMatrixCSR) where {T}
        return tmul2!(y, A, x)
    end

    @eval function  mul!(y::AbstractVector, x::T, A::SparseMatrixCSR) where {T}
        return tmul2!(y, A, x)
    end

    @eval function  *(A::SparseMatrixCSR, x::T) where {T}
        return tmul2(A, x)
    end

    @eval function  *(x::T, A::SparseMatrixCSR) where {T}
        return tmul2(A, x)
    end

end