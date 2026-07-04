# Aqui se aplica la excitacion circular horizontal directamente desde la velocidad en el
# integrador.
using StaticArrays
using LinearAlgebra

mutable struct Particle{N, T}
    r::SVector{N, T}      # Posición
    v::SVector{N, T}      # Velocidad
    a::SVector{N, T}      # Aceleración
    mass::T               # Masa
    radius::T             # Radio
end

function velocidad_marco_inercial(t)
    # Parametros de la excitacion orbital
    amplitud = 5.0                 # Radio de la excitacion
    frecuencia = 1.0                # Frecuencia de la excitacion (Hz)
    omega = 2.0 * pi * frecuencia
    tau = 0.01                        # Tiempo de rampa (s)

    #Factor exponencial para la envolvente de amplitud
    factor_exp = exp(-t / tau)

    # Amplitud y su derivada exacta
    A_t = amplitud * (1.0 - factor_exp * (1.0 + t / tau))
    v_A = amplitud * (t / tau^2) * factor_exp

    # Velocidad exacta 
    v_x = -(v_A * cos(omega * t) - A_t * omega * sin(omega * t))
    v_y = -(v_A * sin(omega * t) + A_t * omega * cos(omega * t))

    return @SVector [v_x, v_y]
end

function velocidad_marco_inercial_simple(t)
    # Parametros de la excitacion orbital
    amplitud = 3.0                 # Radio de la excitacion
    frecuencia = 1.0                # Frecuencia de la excitacion (Hz)
    omega = 2.0 * pi * frecuencia

    # Velocidad exacta 
    v_x = -amplitud * omega * sin(omega * t)
    v_y = amplitud * omega * cos(omega * t)

    return @SVector [v_x, v_y]
end

function velocity_verlet_step!(particles::Vector{Particle{N, T}}, dt::T, tiempo::T, calc_forces!) where {N, T}
    v_marco_actual = velocidad_marco_inercial(tiempo)
    v_marco_siguiente = velocidad_marco_inercial(tiempo + dt)
    delta_v_inercial = v_marco_siguiente - v_marco_actual

    for p in particles
        # r y v evolucionan de forma estandar
        p.r = p.r + p.v * dt + 0.5 * p.a * dt^2
        p.v = p.v + 0.5 * p.a * dt
    end

    calc_forces!(particles, tiempo + dt)

    for p in particles
        p.v = p.v + 0.5 * p.a * dt

        # Se añade la contribución de la aceleración inercial
        p.v += delta_v_inercial
    end
end

function fuerza_total!(particles::Vector{Particle{N, T}}, tiempo::T) where {N, T}
    # Reiniciar aceleraciones
    for p in particles
        p.a = @SVector zeros(T, N)
    end

    # Aplicar fuerzas
    #excitacion_orbital_rampa!(particles, tiempo)
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

    archivo_salida = "resultados/giro_corregido_swirling.xyz"
    open(archivo_salida, "w") do io end # Limpiar archivo si existe

    println("Iniciando movimiento swirled...")

    for paso in 1:pasos
        velocity_verlet_step!(sistema, dt, dt * paso, fuerza_total!)

        if paso % frecuencia_guardado == 0
            guardar_frame_xyz!(archivo_salida, sistema, paso * dt, radio_contenedor)
        end
    end

    println("¡Simulación terminada! Archivo generado: $archivo_salida")
end

simular_sistema()