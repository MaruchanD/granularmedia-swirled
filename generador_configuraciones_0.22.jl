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

function generar_configuracion(particles::Vector{Particle{N, T}}, R_contenedor::Float64; max_iters::Int=20000, dt::Float64=0.1, tol::Float64=1e-6) where {N, T}

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
    R_eff = R_contenedor - maximum(p.radius for p in particles)  # Radio efectivo del recipiente considerando el radio máximo de las partículas
    
    for i in 1:num_p
        p_i = particles[i]                  # Primero se elige a la particula i-esima
        R_eff = R_contenedor - p_i.radius   # Se calcula el radio efectivo del recipiente, considerando el radio de la particula
        r = R_eff * sqrt(rand())            # Se calcula la distancia radial de la particula al centro del recipiente, considerando una distribucion uniforme en el area del circulo
        theta = 2 * pi * rand()             # Se calcula el angulo polar de la particula
        p_i.r = @SVector[r * cos(theta), r * sin(theta)]    # Se asigna la posicion de la particula en coordenadas cartesianas
    end
    println("Las posiciones de las particulas fueron asignadas.")

    # 2. Bucle de Relajacion
    diametro = 2.0 * particles[1].radius
    diam_cuadrado = diametro^2
    desplazamientos = zeros(Float64, num_p, 2)

    for iter in 1:max_iters
        max_superposicion = 0.0
        fill!(desplazamientos, 0.0) # Reiniciar desplazamientos
        
        # Interacción Partícula-Partícula
        for i in 1:(num_p-1)
            for j in (i+1):num_p
                # Primero se calcula la distancia entre las particulas i y j, considerando sus posiciones en coordenadas cartesianas
                dx = particles[i].r[1] - particles[j].r[1]
                dy = particles[i].r[2] - particles[j].r[2]
                dist_cuadrado = dx^2 + dy^2
                
                # Se evalua si las particulas se superponen, considerando el diametro de las particulas. 
                # Si es asi, se calcula la fuerza repulsiva entre ellas y se actualizan los desplazamientos de cada particula
                if dist_cuadrado < diam_cuadrado && dist_cuadrado > 0.0
                    dist = sqrt(dist_cuadrado)
                    superposicion = diametro - dist
                    max_superposicion = max(max_superposicion, superposicion)
                    
                    # Fuerza repulsiva en dirección normal
                    fx = (dx / dist) * superposicion * 0.5
                    fy = (dy / dist) * superposicion * 0.5
                    
                    desplazamientos[i, 1] += fx
                    desplazamientos[i, 2] += fy
                    desplazamientos[j, 1] -= fx
                    desplazamientos[j, 2] -= fy
                end
            end
        end
        
        # Interacción Partícula-Pared
        for i in 1:num_p
            r_sq = particles[i].r[1]^2 + particles[i].r[2]^2
            if r_sq > R_eff^2
                dist = sqrt(r_sq)
                superposicion = dist - R_eff
                max_superposicion = max(max_superposicion, superposicion)

                # Fuerza repulsiva en dirección normal hacia el centro del recipiente
                fx = -(particles[i].r[1] / dist) * superposicion
                fy = -(particles[i].r[2] / dist) * superposicion

                desplazamientos[i, 1] += fx
                desplazamientos[i, 2] += fy
            end
        end
    
        # Actualizar posiciones
        for i in 1:num_p
            particles[i].r += @SVector[desplazamientos[i, 1], desplazamientos[i, 2]] * dt
        end

        # Criterio de convergencia: si la superposición máxima es imperceptible
        if max_superposicion < tol
            println("Convergencia alcanzada en la iteración ", iter)
            break
        end
    end
    #@warn "El sistema no convergió completamente tras $max_iters iteraciones. Superposición máxima residual: $max_superposicion"
end

# ==========================================
# Ejecución y Visualización
# ==========================================

# Hay que hacer modificaciones a la funcion que genera la configuracion para que funcione en 3D. Por ahora, solo funciona en 2D.
# El codigo de visualizacion en 3D es mas complicado, por lo que se dejara para una version futura. Por ahora, solo se visualizara en 2D.
# El numero de particulas junto con el radio del recipiente y el radio de las particulas se pueden modificar para probar diferentes configuraciones.#
# No obstante, hay que tomar en cuenta de no exceder el limite teorico de empaquetamiento en 2D (0.84) para que el sistema pueda converger. En caso de exceder dicho limite, se mostrara un mensaje de advertencia.
# y el programa podria dejar de funcionar.
Dimension = 2
N_particulas = 80 
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