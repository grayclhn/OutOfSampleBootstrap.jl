scalar(a::Number) = a
scalar(A) = A[1]

## Taken from version 0.4
cld{T<:Integer }(x::T, y::T) = div(x,y)+(!signbit(x$y)&(rem(x,y)!=0))
