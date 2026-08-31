# Se utiliza la libreria Plots para visualizar graficamente las posiciones de las particulas del sistema.
using StaticArrays
using LinearAlgebra
using Plots

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

function generar_sistema(N::Int, D::Int; m::Float64=1.0, R::Float64=1.0)
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

# ==========================================
# Ejecución y Visualización
# ==========================================

Dimension = 2
N_particulas = 10
Radio_recipiente = 10.0
Radio_particula = 1.0
Masa_particula = 1.0

# Generar posiciones
sistema = generar_sistema(N_particulas, Dimension, R = Radio_particula)
generar_configuracion(sistema, Radio_recipiente)

# Graficar el resultado
gr() # Usar el backend estándar
plt = plot(aspect_ratio=:equal, legend=false, title="Estado Inicial", grid=true)

# Dibujar el recipiente (aproximado con muchos puntos)
theta = range(0, 2pi, length=200)
plot!(plt, Radio_recipiente .* cos.(theta), Radio_recipiente .* sin.(theta), color=:black, linewidth=2)

# Dibujar las partículas
for i in 1:length(sistema)
    x, y = sistema[i].r
    # plot! en Julia requiere pasar los puntos del contorno para dibujar círculos eficientemente
    plot!(plt, x .+ Radio_particula .* cos.(theta), y .+ Radio_particula .* sin.(theta), fill=(0, :blue), fillalpha=0.6, linecolor=:black)
end

display(plt)