using StaticArrays

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
    # Paso 1 y 2: Actualizar posiciones y velocidades intermedias (medio paso)
    for p in particles
        p.r = p.r + p.v * dt + 0.5 * p.a * dt^2
        p.v = p.v + 0.5 * p.a * dt
    end

    # Paso 3: Calcular las nuevas aceleraciones mutando p.a
    calc_forces!(particles)

    # Paso 4: Completar el paso de velocidad con la nueva aceleración
    for p in particles
        p.v = p.v + 0.5 * p.a * dt
    end
end

# ==========================================
# 3. FUNCIÓN DE FUERZAS TEMPORAL (GRAVEDAD)
# ==========================================
function aplicar_gravedad!(particles::Vector{Particle{N, T}}) where {N, T}
    # Definimos la aceleración de la gravedad según la dimensión
    if N == 2
        g = @SVector[0.0, -9.81]
    else
        g = @SVector[0.0, 0.0, -9.81]
    end

    # Asignamos la gravedad directamente a la aceleración de cada partícula
    for p in particles
        p.a = g
    end
end

# ==========================================
# 4. EJECUCIÓN Y MONITOREO DE LA PRUEBA (1)
# ==========================================
function correr_simulacion()
    # Inicialización de las 3 partículas en 2D usando Float64
    mis_particulas = [
        Particle{2, Float64}(@SVector([0.0, 0.0]), @SVector([1.0, 0.0]), @SVector([0.0, 0.0]), 1.0, 0.5),
        Particle{2, Float64}(@SVector([1.0, 2.0]), @SVector([0.0, 0.0]), @SVector([0.0, 0.0]), 1.5, 0.3),
        Particle{2, Float64}(@SVector([-1.0, 5.0]), @SVector([-0.5, 1.0]), @SVector([0.0, 0.0]), 0.8, 0.4)
    ]

    dt = 0.01  # Paso de tiempo de 10 milisegundos
    pasos = 3  # Ejecutaremos solo 3 pasos para analizar los datos

    println("=== ESTADO INICIAL ===")
    for (i, p) in enumerate(mis_particulas)
        println("Partícula $i -> r: $(p.r), v: $(p.v), a: $(p.a)")
    end

    println("\n=== INICIANDO INTEGRACIÓN ===")
    for paso in 1:pasos
        velocity_verlet_step!(mis_particulas, dt, aplicar_gravedad!)
        
        println("\n--- Paso $paso (t = $(paso * dt) s) ---")
        for (i, p) in enumerate(mis_particulas)
            println("Partícula $i -> r: $(p.r), v: $(p.v), a: $(p.a)")
        end
    end
end
# ==========================================
# 5. EXPORTACIÓN DE DATOS A FORMATO XYZ
# ==========================================
function guardar_frame_xyz!(archivo::String, particles::Vector{Particle{N, T}}, tiempo::Float64) where {N, T}
    # Usamos modo "a" (append) para añadir fotogramas sin borrar los anteriores.
    # El bloque 'do' garantiza que el archivo se cierre correctamente al terminar.
    open(archivo, "a") do io
        # 1. Escribir el número de partículas
        println(io, length(particles))
        
        # 2. Escribir la línea de propiedades y el tiempo actual de la simulación
        println(io, "Properties=species:S:1:pos:R:3:radius:R:1 Time=$tiempo")
        
        # 3. Iterar sobre el arreglo y escribir la data in-place
        for p in particles
            x = p.r[1]
            y = p.r[2]
            # Si el sistema es 2D, forzamos z = 0.0. Si es 3D, tomamos p.r[3]
            z = N == 2 ? 0.0 : p.r[3]
            
            # Escribir la línea: Tipo de partícula ("Granito"), x, y, z, radio
            println(io, "Granito $x $y $z $(p.radius)")
        end
    end
end
# ==========================================
# 6. FUNCIÓN PARA CORRER LA SIMULACIÓN Y EXPORTAR (2)
# ==========================================
function correr_simulacion_con_exportacion()
    # ... (Inicialización de mis_particulas y variables) ...

    mis_particulas = [
        Particle{2, Float64}(@SVector([0.0, 0.0]), @SVector([1.0, 0.0]), @SVector([0.0, 0.0]), 1.0, 0.5),
        Particle{2, Float64}(@SVector([1.0, 2.0]), @SVector([0.0, 0.0]), @SVector([0.0, 0.0]), 1.5, 0.3),
        Particle{2, Float64}(@SVector([-1.0, 5.0]), @SVector([-0.5, 1.0]), @SVector([0.0, 0.0]), 0.8, 0.4)
    ]

    dt = 0.01
    pasos = 100
    archivo_salida = "simulacion_granular.xyz"

    # Limpiar el archivo si ya existe (modo "w" de write)
    open(archivo_salida, "w") do io end 

    println("Iniciando simulación y exportando a $archivo_salida ...")
    
    for paso in 1:pasos
        # 1. Integración matemática
        velocity_verlet_step!(mis_particulas, dt, aplicar_gravedad!)
        
        # 2. Exportar datos (Guardar fotograma)
        tiempo_actual = paso * dt
        guardar_frame_xyz!(archivo_salida, mis_particulas, tiempo_actual)
    end
    
    println("¡Simulación terminada!")
end

# Ejecutar la función principal
# correr_simulacion()
correr_simulacion_con_exportacion()