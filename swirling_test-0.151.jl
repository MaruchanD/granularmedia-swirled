# Las versiones 0.15* agregan modificaciones a las funciones del calculo de la fuerza para
# calcular la energia potencial de la interaccion particula-pared y particula-particula, al igual que la energia cinetica de las particulas.
# Adicionalmente, tambien se buscara agregar la disipacion de energia por friccion. 
# Se agrega ahora las funciones que generan tanto al sistema como una configuracion aleatoria para seguir probando la estabilidad.
# Tambien se hicieron modificaciones mayores para mejorar la modularidad del codigo respecto a versiones anteriores.
# Al igual que se embebio el calculo de la energia potencial en las funciones que calculan los contactos particula-particula y
# particula-pared, para que se pueda calcular la energia mecanica total del sistema de particulas en el laboratorio.

using StaticArrays
using LinearAlgebra

# -- Estructura de datos para representar una particula -- #
mutable struct Particle{N, T}
    # -- Grados de liberta traslacionales -- #
    r::SVector{N, T}      # Posición
    v::SVector{N, T}      # Velocidad
    a::SVector{N, T}      # Aceleración
    mass::T               # Masa
    radius::T             # Radio
end

# Funciones que generan las configuraciones de particulas.

function generar_sistema(N::Int, D::Int; m::Float64=1.0, R::Float64=1.0)
    sistema = [Particle(zeros(SVector{D, Float64}), zeros(SVector{D, Float64}), zeros(SVector{D, Float64}), m, R) for i in 1:N]
    println("El sistema de particulas de dimension ", D," con ", N, " particulas fue creado.")
    return sistema
end

function generar_configuracion!(particles::Vector{Particle{N, T}}, R_contenedor::Float64; max_iters::Int=20000, dt::Float64=0.1, tol::Float64=1e-6) where {N, T}

    num_p = length(particles)

    # Validación física del empaquetamiento en 2D (Límite en 2D es aprox 0.84 - Random Close Packing)
    area_particulas = 0.0
    for i in 1:num_p
        area_particulas += pi * particles[i].radius^2
    end
    area_recipiente = pi * R_contenedor^2
    fraccion_empaquetamiento = area_particulas / area_recipiente
    
    println("Fracción de empaquetamiento requerida: ", round(fraccion_empaquetamiento, digits=3))
    if fraccion_empaquetamiento > 0.82
        @warn "La fracción de empaquetamiento está muy cerca o supera el límite teórico de Random Close Packing en 2D. Es posible que el sistema no logre converger."
    end

    # 1. Inserción inicial aleatoria (se permiten superposiciones)  
    R_eff = R_contenedor - maximum(p.radius for p in particles)  # Radio efectivo del recipiente considerando el radio máximo de las partículas
    
    for i in 1:num_p
        p_i = particles[i]                  # Primero se elige a la particula i-esima
        R_eff = R_contenedor - p_i.radius   # Se calcula el radio efectivo del recipiente, considerando el radio de la particula
        r = R_eff * sqrt(rand())            # Se calcula la distancia radial de la particula al centro del recipiente, considerando una distribucion uniforme en el area del circulo
        theta = 2 * pi * rand()             # Se calcula el angulo polar de la particula
        p_i.r = @SVector[r * cos(theta), r * sin(theta)]    # Se asigna la posicion de la particula en coordenadas cartesianas
    end
    println("Las posiciones de las particulas fueron asignadas.")

    # 2. Bucle de Relajacion
    diametro = 2.0 * particles[1].radius
    diam_cuadrado = diametro^2
    desplazamientos = zeros(Float64, num_p, 2)

    for iter in 1:max_iters
        max_superposicion = 0.0
        fill!(desplazamientos, 0.0) # Reiniciar desplazamientos
        
        # Interacción Partícula-Partícula
        for i in 1:(num_p-1)
            for j in (i+1):num_p
                # Primero se calcula la distancia entre las particulas i y j, considerando sus posiciones en coordenadas cartesianas
                dx = particles[i].r[1] - particles[j].r[1]
                dy = particles[i].r[2] - particles[j].r[2]
                dist_cuadrado = dx^2 + dy^2
                
                # Se evalua si las particulas se superponen, considerando el diametro de las particulas. 
                # Si es asi, se calcula la fuerza repulsiva entre ellas y se actualizan los desplazamientos de cada particula
                if dist_cuadrado < diam_cuadrado && dist_cuadrado > 0.0
                    dist = sqrt(dist_cuadrado)
                    superposicion = diametro - dist
                    max_superposicion = max(max_superposicion, superposicion)
                    
                    # Fuerza repulsiva en dirección normal
                    fx = (dx / dist) * superposicion * 0.5
                    fy = (dy / dist) * superposicion * 0.5
                    
                    desplazamientos[i, 1] += fx
                    desplazamientos[i, 2] += fy
                    desplazamientos[j, 1] -= fx
                    desplazamientos[j, 2] -= fy
                end
            end
        end
        
        # Interacción Partícula-Pared
        for i in 1:num_p
            r_sq = particles[i].r[1]^2 + particles[i].r[2]^2
            if r_sq > R_eff^2
                dist = sqrt(r_sq)
                superposicion = dist - R_eff
                max_superposicion = max(max_superposicion, superposicion)

                # Fuerza repulsiva en dirección normal hacia el centro del recipiente
                fx = -(particles[i].r[1] / dist) * superposicion
                fy = -(particles[i].r[2] / dist) * superposicion

                desplazamientos[i, 1] += fx
                desplazamientos[i, 2] += fy
            end
        end
    
        # Actualizar posiciones
        for i in 1:num_p
            particles[i].r += @SVector[desplazamientos[i, 1], desplazamientos[i, 2]] * dt
        end

        # Criterio de convergencia: si la superposición máxima es imperceptible
        if max_superposicion < tol
            println("Convergencia alcanzada en la iteración ", iter)
            break
        end
    end
    #@warn "El sistema no convergió completamente tras $max_iters iteraciones. Superposición máxima residual: $max_superposicion"
end

# -- Funciones que generan la dinamica del sistema de particulas -- #

# -- Funcion para realizar un paso de integracion de Velocity Verlet -- #
# Se le pasa el arreglo que contiene las particulas, el paso de tiempo, la funcion que calcula las fuerzas, el tiempo actual,
# la tupla que contiene todos los parametros de la simulacion y el radio del recipiente. Se actualizan las posiciones, velocidades y aceleraciones de cada particula.
function velocity_verlet_step!(particles::Vector{Particle{N, T}}, dt::T, calc_forces!, tiempo::T, parametros_Generales::NamedTuple, energia_Mecanica::Vector{T}, radio_Recipiente::T) where {N, T}
    for p in particles
        # -- Calculo de variables de traslacion --#
        p.r = p.r + p.v * dt + 0.5 * p.a * dt^2
        p.v = p.v + 0.5 * p.a * dt
    end
    
    calc_forces!(particles, tiempo + dt, parametros_Generales, energia_Mecanica, radio_Recipiente)
    
    for p in particles
        #-- Calculo de variables de traslacion --#
        p.v = p.v + 0.5 * p.a * dt
    end
end

# -- Funcion para aplicar la excitacion orbital con rampa a las particulas del sistema -- #
function excitacion_orbital_rampa!(particles::Vector{Particle{N,T}}, tiempo::T, (; amplitud, frecuencia, tau)) where {N, T}                         
    
    # Velocidad angular de la excitacion
    omega = 2.0 * pi * frecuencia

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

# -- Funcion para aplicar las fuerzas de contacto entre las particulas y el contenedor circular -- #
function contenedor_circular!(particles::Vector{Particle{N, T}}, radio_Recipiente::T, (; k_n_wall, gamma_n_wall), energia_Mecanica::Vector{T}) where {N, T}
    # Variable para acumular la energia potencial particula-pared
    energia_Potencial_particula_pared = 0.0
    
    # El modelo de fuerzas usado es linear spring-dashpot (resorte + amortiguamiento)
    for p in particles
        d = norm(p.r)                               # Distancia de la particula al centro del recipiente
        delta = d + p.radius - radio_Recipiente     # Solapamiento de la particula con el recipiente

        # Aplicacion de las fuerzas de contacto normales entre particulas y contenedor
        if delta > 0.0
            n_wall = -p.r / d               # Vector normal entre particula y contenedor
            v_n = dot(p.v, n_wall)          # Componente normal de la velocidad de la particula

            # Fuerza de contacto normal con la pared
            F_n_mag = max(k_n_wall * delta - gamma_n_wall * v_n, 0.0)
            F_n_vec = F_n_mag * n_wall

            # Aplicacion de la fuerza normal
            p.a += F_n_vec / p.mass

            # Calculo de la energia potencial de la interaccion particula-pared
            energia_Potencial_particula_pared += 0.5 * k_n_wall * delta^2
        end
    end
    # Registro de la energia potencial particula-pared un instante de tiempo.
    energia_Mecanica[2] = energia_Potencial_particula_pared
end

# -- Funcion para aplicar las fuerzas de contacto entre las particulas -- #
function contacto_particulas!(particles::Vector{Particle{N, T}}, (; k_n, gamma_n), energia_Mecanica::Vector{T}) where {N, T}
    # Variable para acumular la energia potencial particula-particula
    energia_Potencial_particula_particula = 0.0

    # Numero de particulas en el sistema
    num_p = length(particles)

    # Bucle doble para calcular las fuerzas de contacto entre todas las particulas del sistema
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
                v_n_ij = dot(v_ij, n_ij)            # Componente normal de la velocidad relativa

                # Calculo de la fuerzas de contacto normales entre particulas
                # El modelo de fuerzas es linear spring-dashpot
                F_n_mag = max(k_n * delta - gamma_n * v_n_ij, 0.0)
                F_n_vec = F_n_mag * n_ij

                # Aplicacion de la fuerza a las particulas i y j via 3era ley de Newton.
                particles[i].a += F_n_vec / p_i.mass
                particles[j].a -= F_n_vec / p_j.mass

                # Calculo de la energia potencial de la interaccion particula-particula
                energia_Potencial_particula_particula += 0.5 * k_n * delta^2
            end
        end
    end
    # Registro de la energia potencial particula-particula un instante de tiempo.
    energia_Mecanica[3] = energia_Potencial_particula_particula
end

# -- Funcion que calcula todas las fuerzas que actuan sobre las particulas del sistema -- #
function fuerza_total!(particles::Vector{Particle{N, T}}, tiempo::T, parametros_Generales::NamedTuple, energia_Mecanica::Vector{T}, radio_Recipiente::T) where {N, T}
    # Reiniciar aceleraciones
    for p in particles
        p.a = @SVector zeros(T, N)
    end

    # -- Aceleracion inercial asociado al movimiento de swirling del recipiente -- #
    excitacion_orbital_rampa!(particles, tiempo, parametros_Generales)

    # -- Fuerzas debidas al contacto con el recipiente circular -- #
    contenedor_circular!(particles, radio_Recipiente, parametros_Generales, energia_Mecanica)

    # -- Fuerzas debidas al contacto entre particulas -- #
    contacto_particulas!(particles, parametros_Generales, energia_Mecanica)
end

# -- Funcion para calcular la velocidad del marco de referencia del recipiente -- #
function velocidad_marco(tiempo::T, (; amplitud, frecuencia, tau)) where {T}
    
    # Velocidad angular de la excitacion
    omega = 2.0 * pi * frecuencia

    # Amplitud de la excitacion y sus derivadas
    factor_exp = exp(-tiempo/tau) # Factor exponencial
    
    A_t = amplitud * (1. - factor_exp * (1.0 + tiempo/tau))
    v_A = amplitud * (tiempo / tau^2) * factor_exp

    # Componentes de la velocidad del marco de referencia del recipiente.
    v_fx = v_A * cos(omega * tiempo) - A_t * omega * sin(omega * tiempo)
    v_fy = v_A * sin(omega * tiempo) + A_t * omega * cos(omega * tiempo)

    return @SVector [v_fx, v_fy]
end

# -- Funcion para calcular la energia cinetica total del sistema de particulas -- #
# Esta funcion modifica el arreglo de energia_Mecanica, donde el primer elemento es la energia cinetica total del sistema de particulas en el laboratorio.
function calcular_energia_cinetica_laboratorio!(particles::Vector{Particle{N, T}}, tiempo::T, parametros_Generales::NamedTuple, energia_Mecanica::Vector{T}) where {N, T}
    # Primero se calcula la velocidad del marco de referencia del recipiente, y luego se calcula la energia cinetica total del sistema de particulas en el laboratorio.
    V_f = velocidad_marco(tiempo, parametros_Generales)
    E_k = 0.0
    
    for p in particles
        v_abs = p.v + V_f                       # Velocidad absoluta de la particula en el laboratorio
        E_k += 0.5 * p.mass * dot(v_abs, v_abs) # Energia cinetica de la particula en el laboratorio
    end
    
    energia_Mecanica[1] = E_k
end

# -- Funciones para exportar datos de la simulacion -- #

# Exportacion a OVITO
function guardar_frame_xyz(archivo::String, particles::Vector{Particle{N, T}}, tiempo::Float64, radio_contenedor::T) where {N, T}
    open(archivo, "a") do io
        println(io, length(particles)+1)  # Número de partículas + contenedor
        println(io, "Properties=species:S:1:pos:R:3:radius:R:1 Time=$tiempo")
        for p in particles
            println(io, "Granulo $(p.r[1]) $(p.r[2]) 0.0 $(p.radius)")
        end
        println(io, "Contenedor 0.0 0.0 0.0 $radio_contenedor")
    end
end

# -- Funcion para exportar datos extraidos de la simulacion -- #
function guardar_datos(archivo::String, tiempo::Float64, energia_Mecanica::Vector{T}) where {T}
    open(archivo, "a") do io
        println(io, "$tiempo,$(energia_Mecanica[1]),$(energia_Mecanica[2]),$(energia_Mecanica[3]),$(energia_Mecanica[1] + energia_Mecanica[2] + energia_Mecanica[3])")
    end
end

# -- Funcion principal para simular el sistema de particulas -- #
function simular_sistema()
    # -- Parametros del sistema -- #
    dimension_Sistema = 2
    numero_Particulas = 100
    radio_Recipiente = 6.0          # En centimetros
    radio_Particula = 0.5           # En centimetros
    masa_Particula = 10.0            # En gramos

    # -- Parametros para la dinamica del sistema -- #
    # LAS COMAS SON IMPORTANTES PORQUE SON ELEMENTOS DE UNA TUPLA, SI NO SE PONEN LAS COMAS, EL COMPILADOR NO LOS RECONOCE COMO ELEMENTOS DE LA TUPLA.
    parametros_Generales = (
        # -- Parametros de la excitacion orbital -- #
        amplitud = 1.5,     # Radio de la excitacion (1.0 a 5.0 cm)
        frecuencia = 1.0,   # Frecuencia de la excitacion (0.1 - 5.0 Hz)
        tau = 1.0,          # Tiempo de rampa (s)

        # -- Parametros de colisiones particula-recipiente -- #
        
        # -- Normales -- #
        k_n_wall = 1.0e5, #
        gamma_n_wall = 1.0e3, #

        # -- Tangenciales -- #
        k_t_wall = 1.0e5, #
        gamma_t_wall = 1.0e3, #

        # -- Parametros de colisiones particula-particula -- #
        
        # -- Normales -- #
        k_n = 1.0e5, #
        gamma_n = 1.0e3, #

        # -- Tangenciales -- #
        k_t = 1.0e5, #
        gamma_t = 1.0e3, #
    )

    # -- Inicializacion del sistema de particulas -- # 
    sistema = generar_sistema(numero_Particulas, dimension_Sistema, R = radio_Particula, m = masa_Particula)
    generar_configuracion!(sistema, radio_Recipiente)

    # Sistema de particulas de prueba para verificar funciones o estabilidad de la simulacion
    #sistema = [Particle(@SVector[1.0, 0.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5), Particle(@SVector[-1.0, 0.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5)]

    # -- Parametros de la ejecucion -- #
    dt = 1.0e-4                 # En segundos (s)
    tiempo_total = 200.0         # En segundos (s)
    pasos = round(Int, tiempo_total / dt)
    frecuencia_guardado = 100   # Cada 0.01 segundos (s) se guardan los datos de la simulacion

    # -- Archivos de salida -- #
    # Nombre de los archivos de salida para los datos de la simulacion y para los datos extraidos de la simulacion
    # MALDITA SEA, ACUERDATE DE CAMBIAR EL NOMBRE DEL ARCHIVO DE SALIDA, NO TE VUELVAS A EQUIVOCAR
    archivo_salida_1 = "resultados/giro_swirling_151.xyz"
    archivo_salida_2 = "datos/datos_swirling_151.csv"
    # Apertura y limpieza de los archivos de salida
    open(archivo_salida_1, "w") do io end # Limpiar archivo si existe
    open(archivo_salida_2, "w") do io end # Limpiar archivo si existe

    # -- Arreglo para registrar la energia del sistema -- # 
    energia_Mecanica = [
        0.0,    # Energia cinetica total 
        0.0,    # Energia potencial total particula-pared
        0.0     # Energia potencial total particula-particula
        ]

    println("Iniciando movimiento swirled...")
    
    # -- Bucle principal de la simulacion -- #
    for paso in 1:pasos
        # -- Integracion de las ecuaciones de movimiento mediante el metodo de Velocity Verlet -- #
        velocity_verlet_step!(sistema, dt, fuerza_total!, dt * paso, parametros_Generales, energia_Mecanica, radio_Recipiente)
        
        # -- Guardado de los datos de la simulacion y de los datos extraidos de la simulacion -- #
        if paso % frecuencia_guardado == 0
            guardar_frame_xyz(archivo_salida_1, sistema, paso * dt, radio_Recipiente)
            calcular_energia_cinetica_laboratorio!(sistema, paso * dt, parametros_Generales, energia_Mecanica)
            guardar_datos(archivo_salida_2, paso * dt, energia_Mecanica)
        end
    end
    
    println("¡Simulación terminada!\nArchivos generados: $archivo_salida_1 y $archivo_salida_2")
end

simular_sistema()