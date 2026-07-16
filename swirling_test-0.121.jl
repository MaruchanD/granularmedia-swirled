using StaticArrays
using LinearAlgebra

mutable struct Particle{N, T}
    r::SVector{N, T}      # Posición
    v::SVector{N, T}      # Velocidad
    a::SVector{N, T}      # Aceleración
    mass::T               # Masa
    radius::T             # Radio
end

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
    radio_contenedor = 5.0

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
    sistema = [
        Particle(@SVector[0.0, 0.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5),
        Particle(@SVector[2.0, 1.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5),
        Particle(@SVector[-2.5, -1.5], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5)
    ]
    # Parametros
    radio_contenedor = 5.0

    dt = 1.0e-4
    tiempo_total = 10.0 
    pasos = round(Int, tiempo_total / dt)
    frecuencia_guardado = 100 
    
    archivo_salida = "resultados/giro_swirling_121.xyz"
    open(archivo_salida, "w") do io end # Limpiar archivo si existe

    println("Iniciando movimiento swirled...")
    
    for paso in 1:pasos
        velocity_verlet_step!(sistema, dt, fuerza_total!, dt * paso)
        
        if paso % frecuencia_guardado == 0
            guardar_frame_xyz!(archivo_salida, sistema, paso * dt, radio_contenedor)
        end
    end
    
    println("¡Simulación terminada! Archivo generado: $archivo_salida")
end

simular_sistema()