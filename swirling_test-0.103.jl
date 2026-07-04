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

function excitacion_orbital_compensada!(particles::Vector{Particle{N, T}}, tiempo::T) where {N, T}
    # Parametros de la excitacion orbital
    amplitud = 5.0               # Radio de la excitacion
    frecuencia = 1.0             # Frecuencia de la excitacion (Hz)
    omega = 2.0 * pi * frecuencia
    tau = 0.5                       # Tiempo del pulso (s)

    # Factor de decaimiento exponencial
    exp_t = exp(-tiempo / tau)

    # Correcion del eje X (Eliminar choque numerico)
    # El pulso anula la aceleracion inicial para que a_x(0) = 0.
    # Su integral total a lo largo del tiempo es estrictamente 0
    # para no generar derivas a largo plazo.
    pulso_x = amplitud * omega^2 * exp_t * (1.0 - tiempo / tau)

    aceleracion_x = amplitud * omega^2 * cos(omega * tiempo) * exp_t - pulso_x

    # Correcion del eje Y (Eliminar la deriva fisica)
    # El pulso arranca en cero y su area total bajo la curva
    # resta exactamente la constante amplitud*omega que descentra la orbita.
    pulso_y = amplitud * omega * (tiempo / tau^2) * exp_t

    aceleracion_y = amplitud * omega^2 * sin(omega * tiempo) * exp_t - pulso_y

    #Aplicamos todo esto al sistema
    aceleracion_inercial = @SVector [aceleracion_x, aceleracion_y]

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
    excitacion_orbital_compensada!(particles, tiempo)
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
    tiempo_total = 10.0
    pasos = round(Int, tiempo_total / dt)
    frecuencia_guardado = 100

    archivo_salida = "resultados/giro_corregido_swirling_103.xyz"
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