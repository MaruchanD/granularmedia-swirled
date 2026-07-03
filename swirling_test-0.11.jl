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

function fondo_contenedor!(particles::Vector{Particle{N, T}}, tiempo::T) where {N, T}
    # Parametros para la cinematica y las fuerzas contra el fondo
    A = 3.0
    omega = 1.0
    gamma_c = 0.1
    mu_r = 0.5
    gravedad = 9.81

    # Cinematica del contenedor
    # (El contenedor hace un movimiento de remolino / swirling)
    v_contenedor_x = A * omega * cos(omega * tiempo)
    v_contenedor_y = -A * omega * sin(omega * tiempo)
    v_contenedor = @SVector [v_contenedor_x, v_contenedor_y, 0.0]

    a_contenedor_x = -A * omega^2 * sin(omega * tiempo)
    a_contenedor_y = -A * omega^2 * cos(omega * tiempo)
    a_contenedor = @SVector [a_contenedor_x, a_contenedor_y, 0.0]

    # Interaccion con el fondo (Arrastre viscosa y rodadura)
    for p in particles
        v_relativa_fondo = p.v - v_contenedor
        v = norm(v_relativa_fondo)
        v_normal = v_relativa_fondo / (v + 1e-8)    # Evitar división por cero
        
        # Fuerza viscosa (Amortiguacion con el fondo)
        F_viscosa_fondo = -gamma_c * p.mass * v_relativa_fondo
        
        # Friccion de rodadura (Modelo simplificado)
        Fn_gravedad = p.mass * gravedad
        F_rodadura = -mu_r * Fn_gravedad * v_normal
        
        # Componente inercial por rodadura pura (efecto de la aceleración del contenedor) (2/7 masa * aceleración del contenedor)
        F_inercia_rodadura = (2.0 / 7.0) * p.mass * a_contenedor

        p.a += (F_viscosa_fondo + F_rodadura + F_inercia_rodadura) / p.mass
    end
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

# Ejecucion de prueba
function simular_sistema()
    sistema = [
        Particle(@SVector[0.0, 0.0, 0.0], @SVector[0.0, 0.0, 0.0], @SVector[0.0, 0.0, 0.0], 1.0, 0.5),
        Particle(@SVector[1.0, 1.0, 0.0], @SVector[0.0, 0.0, 0.0], @SVector[0.0, 0.0, 0.0], 1.0, 0.5),
        Particle(@SVector[-1.5, 1.5, 0.0], @SVector[0.0, 0.0, 0.0], @SVector[0.0, 0.0, 0.0], 1.0, 0.5)
    ]
    radio_contenedor = 5.0
    dt = 1.0e-4
    tiempo_total = 10.0 
    pasos = round(Int, tiempo_total / dt)
    frecuencia_guardado = 100 
    
    archivo_salida = "fondo_swirling.xyz"
    open(archivo_salida, "w") do io end # Limpiar archivo si existe

    println("Iniciando movimiento swirled...")
    
    for paso in 1:pasos
        velocity_verlet_step!(sistema, dt, fondo_contenedor!, dt * paso)
        
        if paso % frecuencia_guardado == 0
            guardar_frame_xyz!(archivo_salida, sistema, paso * dt, radio_contenedor)
        end
    end
    
    println("¡Simulación terminada! Archivo generado: $archivo_salida")
end

simular_sistema()