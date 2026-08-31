# Primer paso para construir el generador de configuraciones para las simulaciones
# La funcion generar_sistema(N, m, R) genera un vector que contiene las N particulas de masa m y radio R,
# todas inicializadas a 0 en r, v y a.
# Ademas, se le hace una modificacion a la forma en que se muestran los datos de la estructura Particle
# y se manda a imprimir el resultado para verificar que el sistema fue creado adecuadamente.

using StaticArrays
using LinearAlgebra

mutable struct Particle{N, T}
    r::SVector{N, T}      # Posición
    v::SVector{N, T}      # Velocidad
    a::SVector{N, T}      # Aceleración
    mass::T               # Masa
    radius::T             # Radio
end

function generar_sistema(N::Int, m::Float64, R::Float64)
    sistema = [Particle(@SVector[0.0, 0.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], m, R) for i in 1:N]
    return sistema
end

function Base.show(io::IO, p::Particle)
    print(io, "r = ", p.r, ", v = ", p.v, ", a = ", p.a, ", m = ", p.mass, ", R = ", p.radius)
end
objeto = generar_sistema(10, 1.0, 1.0)
println(objeto)