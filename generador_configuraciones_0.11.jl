# Modificacion a la funcion generar_sistema para que considere las dimensiones de las particulas (2D o 3D).
# Tambien se modifico para que use la funcion zeros junto con la estructura SVector para la creacion del
# arreglo que contiene al sistema de particulas con la estructura Particle.
# Se mantiene la forma de mostrar los datos igual que en generador_configuraciones_0.11.jl
using StaticArrays
using LinearAlgebra

mutable struct Particle{N, T}
    r::SVector{N, T}      # Posición
    v::SVector{N, T}      # Velocidad
    a::SVector{N, T}      # Aceleración
    mass::T               # Masa
    radius::T             # Radio
end

function generar_sistema(N::Int, D::Int, m::Float64, R::Float64)
    sistema = [Particle(zeros(SVector{D, Float64}), zeros(SVector{D, Float64}), zeros(SVector{D, Float64}), m, R) for i in 1:N]
    return sistema
end

function Base.show(io::IO, p::Particle)
    print(io, "r = ", p.r, ", v = ", p.v, ", a = ", p.a, ", m = ", p.mass, ", R = ", p.radius)
end
objeto = generar_sistema(10, 2, 1.0, 1.0)
println(objeto)