# Primera union grande de varias funciones
# Se va a unir la funcion que genera las configuraciones con el resto del codigo que produce la dinamica
# swirling. Esto con el objetivo de empezar a probar configuraciones con mas particulas.
# Se va a implementar preliminarmente con la version que solo considera colisiones frontales.

using StaticArrays
using LinearAlgebra

mutable struct Particle{N, T}
    r::SVector{N, T}      # Posición
    v::SVector{N, T}      # Velocidad
    a::SVector{N, T}      # Aceleración
    mass::T               # Masa
    radius::T             # Radio
end

# Funciones que generan las configuraciones de particulas.

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

#Funciones que generan la dinamica del sistema

function velocity_verlet_step!(particles::Vector{Particle{N, T}}, dt::T, calc_forces!, tiempo::T) where {N, T}
    for p in particles
        p.r = p.r + p.v * dt + 0.5 * p.a * dt^2
        p.v = p.v + 0.5 * p.a * dt
    end
    
    calc_forces!(particles, tiempo + dt)
    
    for p in particles
        p.v = p.v + 0.5 * p.a * dt
    end
end

function excitacion_orbital_rampa!(particles::Vector{Particle{N,T}}, tiempo::T) where {N, T}
    # Parametros para la excitacion
    amplitud = 3.0                      # Radio de la excitacion (1.0 a 5.0 cm)
    frecuencia = 1.0                    # Frecuencia de la excitacion (0.1 - 5.0 Hz)
    tau = 1.0                           # Tiempo de rampa (s)
    omega = 2.0 * pi * frecuencia       # Facilidad para construir las cuentas

    # Envolvente para la amplitud y sus derivadas
    factor_exp = exp(-tiempo/tau)       #Factor exponencial
    
    A_t = amplitud * (1. - factor_exp * (1.0 + tiempo/tau))
    v_A = amplitud * (tiempo / tau^2) * factor_exp
    a_A = (amplitud / tau^2) * factor_exp * (1.0 - tiempo/tau)

    # Aceleracion inercial exacta
    aceleracion_x = -(a_A * cos(omega * tiempo) - 2.0 * v_A * (omega) * sin(omega * tiempo) - A_t * (omega)^2 * cos(omega * tiempo))
    aceleracion_y = -(a_A * sin(omega * tiempo) + 2.0 * v_A * (omega) * cos(omega * tiempo) - A_t * (omega)^2 * sin(omega * tiempo))
    aceleracion_inercial = @SVector [aceleracion_x, aceleracion_y]

    # Aplicar la aceleracion inercial a cada particula
    for p in particles
        p.a += aceleracion_inercial
    end
end

function contenedor_circular!(particles::Vector{Particle{N, T}}) where {N, T}
    # Parametros de la interaccion con el contenedor
    kn_wall = 1.0e5
    gamma_wall = 5.0e3
    radio_contenedor = 10.0

    # Aplicacion de la fuerza de contacto normal con el contenedor.
    # El modelo de fuerzas usado es linear spring-dashpot (resorte + amortiguamiento)
    for p in particles
        d = norm(p.r)
        delta = d + p.radius - radio_contenedor

        if delta > 0.0
            n_wall = -p.r / d
            vn = dot(p.v, n_wall)

            Fn_mag = max(kn_wall * delta - gamma_wall * vn, 0.0)
            p.a += (Fn_mag * n_wall) / p.mass
        end
    end

    # Aplicacion de la fuerza tangencial (friccion) con el contenedor.
end

function contacto_particulas!(particles::Vector{Particle{N, T}}) where {N, T}
    # Parametros de la interaccion entre particulas
    kn = 1.0e5
    gamma_n = 5.0e1
    num_p = length(particles)

    for i in 1:num_p
        p_i = particles[i]
        for j in (i+1):num_p
            p_j = particles[j]

            # Definicion de la distancia entre particulas y la suma de sus radios.
            r_ij = p_i.r - p_j.r
            d = norm(r_ij)
            suma_radios = p_i.radius + p_j.radius

            if d < suma_radios && d > 0.0
                delta = suma_radios - d             # Solapamiento de las particulas
                n_ij = r_ij / d                     # Vector normal entre particulas
                v_ij = p_i.v - p_j.v                # Velocidad relativa entre particulas
                vn = dot(v_ij, n_ij)                # Componente normal de la velocidad relativa

                # Calculo de la fuerzas de contacto normales entre particulas
                # El modelo de fuerzas es linear spring-dashpot
                Fn_mag = max(kn * delta - gamma_n * vn, 0.0)
                Fn_vec = Fn_mag * n_ij
                
                # Aplicacion de la fuerza a las particulas i y j via 3era ley de Newton.
                particles[i].a += Fn_vec / p_i.mass
                particles[j].a -= Fn_vec / p_j.mass

                # Calculo de las fuerzas de contacto tangenciales entre particulas
            end
        end
    end
end

function fuerza_total!(particles::Vector{Particle{N, T}}, tiempo::T) where {N, T}
    # Reiniciar aceleraciones
    for p in particles
        p.a = @SVector zeros(T, N)
    end

    # Aplicar todas las fuerzas
    excitacion_orbital_rampa!(particles, tiempo)
    contenedor_circular!(particles)
    contacto_particulas!(particles)
end

# Exportacion a OVITO
function guardar_frame_xyz!(archivo::String, particles::Vector{Particle{N, T}}, tiempo::Float64, radio_contenedor::T) where {N, T}
    open(archivo, "a") do io
        println(io, length(particles)+1)  # Número de partículas + contenedor
        println(io, "Properties=species:S:1:pos:R:3:radius:R:1 Time=$tiempo")
        for p in particles
            println(io, "Granulo $(p.r[1]) $(p.r[2]) 0.0 $(p.radius)")
        end
        println(io, "Contenedor 0.0 0.0 0.0 $radio_contenedor")
    end
end

function simular_sistema()
    # Parametros del sistema
    dimension_Sistema = 2               # Entero
    n_Particulas = 50                  # Entero
    radio_Recipiente = 10.0             # Real
    radio_Particula = 1.0
    masa_Particula = 1.0

    #Funciones para obtener el sistema
    sistema = generar_sistema(n_Particulas, dimension_Sistema)
    generar_configuracion(sistema, radio_Recipiente)

    # Parametros de la ejecucion
    dt = 1.0e-4
    tiempo_total = 10.0 
    pasos = round(Int, tiempo_total / dt)
    frecuencia_guardado = 100 
    
    archivo_salida = "resultados/giro_swirling_122.xyz" # IMPORTANTE MODIFICAR AL CREAR UN NUEVO ARCHIVO
    open(archivo_salida, "w") do io end # Limpiar archivo si existe

    println("Iniciando movimiento swirled...")
    
    for paso in 1:pasos
        velocity_verlet_step!(sistema, dt, fuerza_total!, dt * paso)
        
        if paso % frecuencia_guardado == 0
            guardar_frame_xyz!(archivo_salida, sistema, paso * dt, radio_Recipiente)
        end
    end
    
    println("¡Simulación terminada! Archivo generado: $archivo_salida")
end

simular_sistema()