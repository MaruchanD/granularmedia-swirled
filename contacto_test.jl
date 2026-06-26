using StaticArrays
using LinearAlgebra

# ==========================================
# 1. DEFINICIÓN DE LA ESTRUCTURA
# ==========================================
mutable struct Particle{N, T}
    r::SVector{N, T}      # Posición (m)
    v::SVector{N, T}      # Velocidad (m/s)
    a::SVector{N, T}      # Aceleración (m/s²)
    mass::T               # Masa (kg)
    radius::T             # Radio (m)
end

# ==========================================
# 2. ALGORITMO DE INTEGRACIÓN (VELOCITY-VERLET)
# ==========================================
function velocity_verlet_step!(particles::Vector{Particle{N, T}}, dt::T, calc_forces!) where {N, T}
    # Paso 1 y 2: Actualizar posiciones y velocidades a medio paso
    for p in particles
        p.r = p.r + p.v * dt + 0.5 * p.a * dt^2
        p.v = p.v + 0.5 * p.a * dt
    end

    # Paso 3: Calcular las nuevas fuerzas/aceleraciones basadas en las nuevas posiciones
    calc_forces!(particles)

    # Paso 4: Completar el paso de velocidad
    for p in particles
        p.v = p.v + 0.5 * p.a * dt
    end
end

# ==========================================
# 3. FUNCIÓN DE FUERZAS DE CONTACTO NORMAL
# ==========================================
function calcular_fuerzas_contacto!(particles::Vector{Particle{N, T}}) where {N, T}
    # En esta prueba no hay gravedad, inicializamos las aceleraciones en cero
    for p in particles
        p.a = zeros(SVector{N, T})
    end

    # Parámetros del modelo mecánico lineal
    kn = 1.0e5  # Rigidez del contacto (N/m)
    gn = 5.0e1  # Coeficiente de amortiguamiento (kg/s) -> Controla la ineslasticidad
    
    num_p = length(particles)

    for i in 1:num_p
        p_i = particles[i]
        for j in (i+1):num_p
            p_j = particles[j]

            # Vector de distancia relativa
            r_ij = p_i.r - p_j.r
            d = norm(r_ij)
            
            suma_radios = p_i.radius + p_j.radius
            
            # Condición de contacto geométrico
            if d < suma_radios && d > 0.0
                delta = suma_radios - d  # Solapamiento
                n_ij = r_ij / d          # Vector normal unitario
                
                # Velocidad relativa entre las dos partículas
                v_ij = p_i.v - p_j.v
                vn = dot(v_ij, n_ij)     # Proyección normal
                
                # Ecuación de fuerza normal: Resorte - Amortiguador
                Fn_mag = kn * delta - gn * vn
                Fn_mag = max(Fn_mag, 0.0) # Evitar fuerzas atractivas artificiales
                
                Fn_vec = Fn_mag * n_ij
                
                # Aplicación de la Tercera Ley de Newton (Acción y Reacción)
                particles[i].a += Fn_vec / p_i.mass
                particles[j].a -= Fn_vec / p_j.mass
            end
        end
    end
end

# ==========================================
# 4. EXPORTACIÓN DE DATOS PARA OVITO
# ==========================================
function guardar_frame_xyz!(archivo::String, particles::Vector{Particle{N, T}}, tiempo::Float64) where {N, T}
    open(archivo, "a") do io
        println(io, length(particles))
        println(io, "Properties=species:S:1:pos:R:3:radius:R:1 Time=$tiempo")
        for p in particles
            x = p.r[1]
            y = p.r[2]
            z = N == 2 ? 0.0 : p.r[3]
            println(io, "Granulo $x $y $z $(p.radius)")
        end
    end
end

# ==========================================
# 5. FUNCIÓN PRINCIPAL DE EJECUCIÓN
# ==========================================
function simular_colision()
    # Inicialización: Dos partículas idénticas dirigiéndose una hacia la otra
    # Partícula 1: Viene desde la izquierda (x = 0.0) hacia la derecha (vx = 2.0)
    # Partícula 2: Viene desde la derecha (x = 1.8) hacia la izquierda (vx = -2.0)
    # Particula 3: Viene desde la derecha (x = 2.5) hacia la izquierda (vx = -1.0)
    # Ambas tienen radio = 0.5, por lo que chocarán cuando la distancia entre centros sea < 1.0
    mis_particulas = [
        Particle{2, Float64}(@SVector([0.0, 0.0]), @SVector([2.0, 0.0]), @SVector([0.0, 0.0]), 5.0, 0.5),
        Particle{2, Float64}(@SVector([1.8, 0.0]), @SVector([-2.0, 0.0]), @SVector([0.0, 0.0]), 1.0, 0.5)
        #,Particle{2, Float64}(@SVector([3.0, 0.0]), @SVector([-3.0, 0.0]), @SVector([0.0, 0.0]), 1.0, 0.5)  # Partícula extra para observar interacciones
    ]

    dt = 1.0e-4          # Paso de tiempo pequeño para mantener la estabilidad del resorte
    tiempo_total = 0.5   # Simular medio segundo en total
    pasos = round(Int, tiempo_total / dt)
    
    # Para no saturar el disco duro, guardamos 1 frame cada 50 pasos de integración
    frecuencia_guardado = 50 
    archivo_salida = "colision_granular.xyz"

    # Resetear el archivo de texto antes de empezar
    open(archivo_salida, "w") do io end 

    println("Ejecutando $pasos pasos de simulación...")
    
    for paso in 1:pasos
        velocity_verlet_step!(mis_particulas, dt, calcular_fuerzas_contacto!)
        
        if paso % frecuencia_guardado == 0
            guardar_frame_xyz!(archivo_salida, mis_particulas, paso * dt)
        end
    end
    
    println("¡Simulación completada con éxito!")
    println("Archivo generado: $archivo_salida")
    println("\n=== ESTADO FINAL DE LAS PARTÍCULAS ===")
    for (i, p) in enumerate(mis_particulas)
        println("Partícula $i -> r: $(p.r), v: $(p.v)")
    end
end

# Ejecutar la simulación
simular_colision()