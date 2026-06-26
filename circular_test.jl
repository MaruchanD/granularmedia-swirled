using StaticArrays
using LinearAlgebra

# ==========================================
# 1. DEFINICIÓN DE LA ESTRUCTURA
# ==========================================
mutable struct Particle{N, T}
    r::SVector{N, T}      # Posición
    v::SVector{N, T}      # Velocidad
    a::SVector{N, T}      # Aceleración
    mass::T               # Masa
    radius::T             # Radio
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
function aplicar_fuerzas_contenedor_circular!(particles::Vector{Particle{N, T}}, radio_contenedor::T) where {N, T}
    kn_wall = 1.0e5
    gn_wall = 5.0e2

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

# Función orquestadora
function calcular_fuerzas_totales!(particles::Vector{Particle{2, Float64}})
    # 1. Resetear aceleración e imponer Gravedad
    g = @SVector[0.0, -9.81]
    for p in particles
        p.a = g
    end

    # 2. Fuerzas del contenedor circular (Radio = 5.0)
    aplicar_fuerzas_contenedor_circular!(particles, 5.0)

    # 3. Fuerzas de contacto inter-particulares
    calcular_fuerzas_contacto!(particles)
end

# ==========================================
# 4. EXPORTACIÓN A OVITO
# ==========================================
function guardar_frame_xyz!(archivo::String, particles::Vector{Particle{N, T}}, tiempo::Float64) where {N, T}
    open(archivo, "a") do io
        println(io, length(particles)+1)  # Número de partículas + contenedor
        println(io, "Properties=species:S:1:pos:R:3:radius:R:1 Time=$tiempo")
        for p in particles
            println(io, "Granulo $(p.r[1]) $(p.r[2]) 0.0 $(p.radius)")
        end
        println(io, "Contenedor 0.0 0.0 0.0 5.0")
    end
end

# ==========================================
# 5. EJECUCIÓN DE LA PRUEBA
# ==========================================
function simular_contenedor_circular()
    # Inicializamos 5 partículas dentro del radio de 5.0 metros
    mis_particulas = [
        Particle{2, Float64}(@SVector([ 0.0,  3.0]), @SVector([0.0, 0.0]), @SVector([0.0, 0.0]), 1.0, 0.5),
        Particle{2, Float64}(@SVector([-1.0,  2.0]), @SVector([0.0, 0.0]), @SVector([0.0, 0.0]), 1.0, 0.5),
        Particle{2, Float64}(@SVector([ 1.0,  1.0]), @SVector([0.0, 0.0]), @SVector([0.0, 0.0]), 1.0, 0.5),
        Particle{2, Float64}(@SVector([-0.5,  4.0]), @SVector([0.0, 0.0]), @SVector([0.0, 0.0]), 1.0, 0.5),
        Particle{2, Float64}(@SVector([ 0.5,  2.5]), @SVector([0.0, 0.0]), @SVector([0.0, 0.0]), 1.0, 0.5)
    ]

    dt = 1.0e-4
    tiempo_total = 10.0 # 10.0 segundos de simulación física
    pasos = round(Int, tiempo_total / dt)
    frecuencia_guardado = 100 
    
    archivo_salida = "contenedor_circular.xyz"
    open(archivo_salida, "w") do io end # Limpiar archivo si existe

    println("Iniciando caída de partículas en el contenedor circular...")
    
    for paso in 1:pasos
        velocity_verlet_step!(mis_particulas, dt, calcular_fuerzas_totales!)
        
        if paso % frecuencia_guardado == 0
            guardar_frame_xyz!(archivo_salida, mis_particulas, paso * dt)
        end
    end
    
    println("¡Simulación terminada! Archivo generado: $archivo_salida")
end

simular_contenedor_circular()