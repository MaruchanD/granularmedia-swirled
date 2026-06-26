using StaticArrays
using LinearAlgebra

# ==========================================
# 1. DEFINICIÓN DE LA ESTRUCTURA
# ==========================================
mutable struct Particle{N, T}
    r::SVector{N, T}
    v::SVector{N, T}
    a::SVector{N, T}
    mass::T
    radius::T
end

# ==========================================
# 2. INTEGRADOR (VELOCITY-VERLET)
# ==========================================
function velocity_verlet_step!(particles::Vector{Particle{N, T}}, dt::T, calc_forces!) where {N, T}
    for p in particles
        p.r = p.r + p.v * dt + 0.5 * p.a * dt^2
        p.v = p.v + 0.5 * p.a * dt
    end
    
    calc_forces!(particles)
    
    for p in particles
        p.v = p.v + 0.5 * p.a * dt
    end
end

# ==========================================
# 3. MÓDULOS DE FUERZAS
# ==========================================
function aplicar_fuerzas_paredes!(particles::Vector{Particle{N, T}}, x_min::T, x_max::T, y_min::T, y_max::T) where {N, T}
    kn_wall = 1.0e5
    gn_wall = 1.0e2  # Aumentamos un poco la disipación para que se asienten más rápido

    for p in particles
        # Pared Izquierda
        if (delta = x_min + p.radius - p.r[1]) > 0.0
            Fn = max(kn_wall * delta - gn_wall * p.v[1], 0.0)
            p.a += SVector{N, T}(Fn / p.mass, 0.0)
        end
        # Pared Derecha
        if (delta = p.r[1] + p.radius - x_max) > 0.0
            Fn = max(kn_wall * delta - gn_wall * (-p.v[1]), 0.0)
            p.a += SVector{N, T}(-Fn / p.mass, 0.0)
        end
        # Suelo
        if (delta = y_min + p.radius - p.r[2]) > 0.0
            Fn = max(kn_wall * delta - gn_wall * p.v[2], 0.0)
            p.a += SVector{N, T}(0.0, Fn / p.mass)
        end
        # Techo
        if (delta = p.r[2] + p.radius - y_max) > 0.0
            Fn = max(kn_wall * delta - gn_wall * (-p.v[2]), 0.0)
            p.a += SVector{N, T}(0.0, -Fn / p.mass)
        end
    end
end

function calcular_fuerzas_contacto!(particles::Vector{Particle{N, T}}) where {N, T}
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

# Función orquestadora de fuerzas
function calcular_fuerzas_totales!(particles::Vector{Particle{2, Float64}})
    # 1. Resetear aceleración e imponer Gravedad
    g = @SVector[0.0, -9.81]
    for p in particles
        p.a = g
    end

    # 2. Fuerzas de las paredes (Silo de ancho 4, alto 10)
    aplicar_fuerzas_paredes!(particles, -2.0, 2.0, 0.0, 10.0)

    # 3. Fuerzas de contacto inter-particulares
    calcular_fuerzas_contacto!(particles)
end

# ==========================================
# 4. EXPORTACIÓN A OVITO
# ==========================================
function guardar_frame_xyz!(archivo::String, particles::Vector{Particle{N, T}}, tiempo::Float64) where {N, T}
    open(archivo, "a") do io
        println(io, length(particles))
        println(io, "Properties=species:S:1:pos:R:3:radius:R:1 Time=$tiempo")
        for p in particles
            println(io, "Granulo $(p.r[1]) $(p.r[2]) 0.0 $(p.radius)")
        end
    end
end

# ==========================================
# 5. EJECUCIÓN DE LA PRUEBA
# ==========================================
function simular_confinamiento()
    # Inicializamos 5 partículas a diferentes alturas para evitar solapamientos iniciales
    mis_particulas = [
        Particle{2, Float64}(@SVector([ 0.0,  2.0]), @SVector([0.0, 0.0]), @SVector([0.0, 0.0]), 1.0, 0.4),
        Particle{2, Float64}(@SVector([-1.0,  4.0]), @SVector([0.0, 0.0]), @SVector([0.0, 0.0]), 1.2, 0.5),
        Particle{2, Float64}(@SVector([ 1.0,  6.0]), @SVector([0.0, 0.0]), @SVector([0.0, 0.0]), 0.8, 0.35),
        Particle{2, Float64}(@SVector([-0.5,  8.0]), @SVector([0.0, 0.0]), @SVector([0.0, 0.0]), 1.5, 0.45),
        Particle{2, Float64}(@SVector([ 0.5, 10.0]), @SVector([0.0, 0.0]), @SVector([0.0, 0.0]), 1.0, 0.4)
    ]

    dt = 1.0e-4
    tiempo_total = 1.5 # 1.5 segundos de simulación física
    pasos = round(Int, tiempo_total / dt)
    frecuencia_guardado = 100 # Guardar frame cada 100 pasos (0.01s en la animación)
    
    archivo_salida = "silo_granular.xyz"
    open(archivo_salida, "w") do io end # Limpiar archivo

    println("Iniciando caída de partículas en el contenedor...")
    
    for paso in 1:pasos
        velocity_verlet_step!(mis_particulas, dt, calcular_fuerzas_totales!)
        
        if paso % frecuencia_guardado == 0
            guardar_frame_xyz!(archivo_salida, mis_particulas, paso * dt)
        end
    end
    
    println("¡Simulación terminada! Abre el archivo '$archivo_salida' en OVITO.")
end

simular_confinamiento()