# Implementacion de una funcion para la excitacion orbital de un contenedor circular,
# con una rampa de amplitud y frecuencia.
# Se construyeron dos funciones, una simple y otra con rampa.
# La simple necesita una funcion que ajuste adecuadamente las velocidades iniciales de las
# particulas, mientras que la de rampa no necesita de esto, ya que la aceleracion inercial
# es exacta y no depende de las velocidades iniciales, ademas de imita el inicio del
# del movimiento de un shaker. No obstante, produce deriva para valores de tau comparables
# con el paso de integracion dt.
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

function excitacion_orbital_simple!(particles::Vector{Particle{N, T}}, tiempo::T) where {N, T}
    # Parametros de la excitacion orbital
    amplitud = 3.0                 # Radio de la excitacion
    frecuencia = 1.0                # Frecuencia de la excitacion (Hz)

    # Calcular el vector de aceleración inercial (Fuerza ficticia / masa)
    aceleracion_inercial_x = amplitud * (2 * pi * frecuencia)^2 * cos(2 * pi * frecuencia * tiempo)
    aceleracion_inercial_y = amplitud * (2 * pi * frecuencia)^2 * sin(2 * pi * frecuencia * tiempo)
    
    aceleracion_inercial = @SVector [aceleracion_inercial_x, aceleracion_inercial_y]

    # Aplicar la aceleración inercial a cada partícula
    for p in particles
        p.a += aceleracion_inercial
    end
end

function excitacion_orbital_rampa!(particles::Vector{Particle{N, T}}, tiempo::T) where {N, T}
    amplitud = 5.0                 # Radio de la excitacion
    frecuencia = 1.0                # Frecuencia de la excitacion (Hz)
    tau = 1.0                        # Tiempo de rampa (s)
    omega = 2.0 * pi * frecuencia

    # Envolvente de amplitud y sus derivadas
    factor_exp = exp(-tiempo / tau)
    A_t = amplitud * (1.0 - factor_exp * (1.0 + tiempo / tau))
    v_A = amplitud * (tiempo / tau^2) * factor_exp
    a_A = (amplitud / tau^2) * factor_exp * (1.0 - tiempo / tau)

    # Aceleracion inercial exacta 
    aceleracion_x = -(a_A * cos(omega * tiempo) - 2.0 * v_A * (omega) * sin(omega * tiempo) - A_t * (omega)^2 * cos(omega * tiempo))
    aceleracion_y = -(a_A * sin(omega * tiempo) + 2.0 * v_A * (omega) * cos(omega * tiempo) - A_t * (omega)^2 * sin(omega * tiempo))
    aceleracion_inercial = @SVector [aceleracion_x, aceleracion_y]

    # Aplicar la aceleración inercial a cada partícula
    for p in particles
        p.a += aceleracion_inercial
    end
end

function fuerza_total!(particles::Vector{Particle{N, T}}, tiempo::T) where {N, T}
    # Reiniciar aceleraciones
    for p in particles
        p.a = @SVector zeros(T, N)
    end

    # Aplicar fuerzas
    excitacion_orbital_rampa!(particles, tiempo)
end

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
        Particle(@SVector[0.0, 0.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5)#,
        #Particle(@SVector[1.0, 0.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5),
        #Particle(@SVector[0.0, 1.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5),
        #Particle(@SVector[-1.0, 0.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5)
    ]
    radio_contenedor = 5.0

    dt = 1.0e-4
    tiempo_total = 20.0
    pasos = round(Int, tiempo_total / dt)
    frecuencia_guardado = 100

    archivo_salida = "resultados/giro_corregido_swirling_101.xyz"
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