# Solo se aplica una fuerza inercial debido al giro del contenedor. 
# Se toman en cuenta interacciones entre particulas y con las paredes del contenedor.
# La funcion que calcula las fuerzas resetea las aceleraciones de las particulas, 
# por lo que no es necesario hacerlo en el integrador. Para evitar acumulacion de fuerzas.
using StaticArrays
using LinearAlgebra

# DEFINICIÓN DE LA ESTRUCTURA DE LA PARTICULA
mutable struct Particle{N, T}
    r::SVector{N, T}      # Posición
    v::SVector{N, T}      # Velocidad
    a::SVector{N, T}      # Aceleración
    mass::T               # Masa
    radius::T             # Radio
end

# INTEGRADOR (VELOCITY-VERLET)
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

# FUERZAS APLICADAS A LAS PARTICULAS
function contenedor_circular!(particles::Vector{Particle{N, T}}) where {N, T}
    kn_wall = 1.0e5
    gn_wall = 0.0 #5.0e3
    radio_contenedor = 5.0

    for p in particles
        d = norm(p.r)
        delta = d + p.radius - radio_contenedor
        
        if delta > 0.0
            n_wall = -p.r / d
            vn = dot(p.v, n_wall)
            
            Fn_mag = max(kn_wall * delta - gn_wall * vn, 0.0)
            p.a += (Fn_mag * n_wall) / p.mass
        end
    end
end

function fuerza_inercial!(particles::Vector{Particle{N, T}}, tiempo::T) where {N, T}
    omega = 2.0 * pi                                      # Velocidad angular del contenedor (rad/s)
    r_giro = 4.0
    for p in particles
        # Fuerza inercial debido al giro del envase
        F_x = p.mass * r_giro * omega^2 * cos(omega * tiempo)
        F_y = p.mass * r_giro * omega^2 * sin(omega * tiempo)
        F_i = @SVector [F_x, F_y]
        # Aceleración centrífuga: a_c = F_c / m = omega^2 * r
        p.a += F_i / p.mass
    end
end

function contacto_particulas!(particles::Vector{Particle{N, T}}) where {N, T}
    kn = 1.0e5
    gn = 5.0e1
    num_p = length(particles)

    for i in 1:num_p
        p_i = particles[i]
        for j in (i+1):num_p
            p_j = particles[j]
            r_ij = p_i.r - p_j.r
            d = norm(r_ij)
            suma_radios = p_i.radius + p_j.radius
            
            if d < suma_radios && d > 0.0
                delta = suma_radios - d
                n_ij = r_ij / d
                v_ij = p_i.v - p_j.v
                vn = dot(v_ij, n_ij)
                
                Fn_mag = max(kn * delta - gn * vn, 0.0)
                Fn_vec = Fn_mag * n_ij
                
                particles[i].a += Fn_vec / p_i.mass
                particles[j].a -= Fn_vec / p_j.mass
            end
        end
    end
end

function fuerza_total!(particles::Vector{Particle{N, T}}, tiempo::T) where {N, T}
    # Resetear aceleraciones 
    # (ESTE PASO ES CRUCIAL PARA EVITAR ACUMULACIÓN DE FUERZAS)
    for p in particles
        p.a = @SVector zeros(T, N)
    end
    contacto_particulas!(particles)                               # Fuerza de contacto entre partículas
    contenedor_circular!(particles)                     # Fuerza ejercida por el contenedor circular
    fuerza_inercial!(particles, tiempo)                         # Fuerza inercial debido al giro del envase
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

# Ejecucion de prueba
function simular_sistema()
    sistema = [
        Particle(@SVector[0.0, 0.0], @SVector[0.0, -8.0 * pi], @SVector[0.0, 0.0], 1.0, 0.5),
        Particle(@SVector[1.0, 1.0], @SVector[0.0, -8.0 * pi], @SVector[0.0, 0.0], 1.0, 0.5),
        Particle(@SVector[-1.5, 1.5], @SVector[0.0, -8.0 * pi], @SVector[0.0, 0.0], 1.0, 0.5)
    ]
    # Parametros
    radio_contenedor = 5.0

    dt = 1.0e-4
    tiempo_total = 5.0 
    pasos = round(Int, tiempo_total / dt)
    frecuencia_guardado = 100 
    
    archivo_salida = "resultados/contenedor_swirling.xyz"
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
