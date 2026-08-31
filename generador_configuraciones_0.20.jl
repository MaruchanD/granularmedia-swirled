# Se agrego la primera version de la funcion generar_configuracion(SISTEMA, R_contenedor). Dicha funcion
# verifica el empaquetamiento fisico de las particulas, al calcular la fraccion de empaquetamiento. Siempre 
# que dicha fraccion sea menor que 0.84 (El limite para un Random Close Packing en 2D), el empaquetamiento
# es posible. 
# Seguidamente, comienza la asignacion aleatoria de posiciones a las particulas del sistema. Como primera
# version, no se esta tomando en cuenta si se superponen o no. Por ultimo, se manda a imprimir el arreglo
# del sistema para verificar que las posiciones fueron asignadas.
using StaticArrays
using LinearAlgebra

mutable struct Particle{N, T}
    r::SVector{N, T}      # Posición
    v::SVector{N, T}      # Velocidad
    a::SVector{N, T}      # Aceleración
    mass::T               # Masa
    radius::T             # Radio
end

function Base.show(io::IO, p::Particle)
    print(io, "r = ", p.r, ", v = ", p.v, ", a = ", p.a, ", m = ", p.mass, ", R = ", p.radius)
end

function generar_sistema(N::Int, D::Int, m::Float64, R::Float64)
    sistema = [Particle(zeros(SVector{D, Float64}), zeros(SVector{D, Float64}), zeros(SVector{D, Float64}), m, R) for i in 1:N]
    println("El sistema de particulas de dimension ", D," con ", N, " particulas fue creado.")
    return sistema
end

function generar_configuracion(particles::Vector{Particle{N, T}}, R_contenedor::Float64) where {N, T}

    num_p = length(particles)

    # Validación física del empaquetamiento en 2D (Límite en 2D es aprox 0.84 - Random Close Packing)
    area_particulas = 0.0
    for i in 1:num_p
        area_particulas += pi * particles[i].radius^2
    end
    area_recipiente = pi * R_contenedor^2
    fraccion_empaquetamiento = area_particulas / area_recipiente
    
    println("Fracción de empaquetamiento requerida: ", round(fraccion_empaquetamiento, digits=3))
    if fraccion_empaquetamiento > 0.82
        @warn "La fracción de empaquetamiento está muy cerca o supera el límite teórico de Random Close Packing en 2D. Es posible que el sistema no logre converger."
    end

    # 1. Inserción inicial aleatoria (se permiten superposiciones)    
    for i in 1:num_p
        p_i = particles[i]
        R_eff = R_contenedor - p_i.radius
        r = R_eff * sqrt(rand())
        theta = 2 * pi * rand()
        p_i.r = @SVector[r * cos(theta), r * sin(theta)]
    end
    println("Las posiciones de las particulas fueron asignadas.")
end

sistema = generar_sistema(10, 2, 1.0, 1.0)
generar_configuracion(sistema, 5)
println(sistema)