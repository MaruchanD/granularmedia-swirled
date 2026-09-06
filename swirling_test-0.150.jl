# Las versiones 0.15* agregan modificaciones a las funciones del calculo de la fuerza para
# calcular la energia potencial de la interaccion particula-pared y particula-particula, al igual que la energia cinetica de las particulas.
# Adicionalmente, tambien se buscara agregar la disipacion de energia por friccion. 
# En este archivo solo integrara la energia cinetica total del sistema, sin considerar la rotacion de las particulas.
# Se lleva a acabo con un sistema de una particula para pruebas simples de estabilidad.

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

# -- Funciones que generan la dinamica del sistema de particulas -- #

# -- Funcion para realizar un paso de integracion de Velocity Verlet -- #

function velocity_verlet_step!(particles::Vector{Particle{N, T}}, dt::T, calc_forces!, tiempo::T, radio_Recipiente::T) where {N, T}
    for p in particles
        # -- Calculo de variables de traslacion --#
        p.r = p.r + p.v * dt + 0.5 * p.a * dt^2
        p.v = p.v + 0.5 * p.a * dt
    end
    
    calc_forces!(particles, tiempo + dt, radio_Recipiente)
    
    for p in particles
        #-- Calculo de variables de traslacion --#
        p.v = p.v + 0.5 * p.a * dt
    end
end

# -- Funcion para aplicar la excitacion orbital con rampa a las particulas del sistema -- #
function excitacion_orbital_rampa!(particles::Vector{Particle{N,T}}, tiempo::T) where {N, T}
    # Parametros para la excitacion
    amplitud = 1.5                      # Radio de la excitacion (1.0 a 5.0 cm)
    frecuencia = 1.0                    # Frecuencia de la excitacion (0.1 - 5.0 Hz)
    tau = 1.0                           # Tiempo de rampa (s)
    omega = 2.0 * pi * frecuencia       # Facilidad para construir las cuentas

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
function contenedor_circular!(particles::Vector{Particle{N, T}}, radio_Recipiente::T) where {N, T}
    # Parametros de la interaccion con el contenedor
    kn_wall = 1.0e5
    gamma_wall = 5.0e3

    # Aplicacion de la fuerza de contacto normal con el contenedor.
    # El modelo de fuerzas usado es linear spring-dashpot (resorte + amortiguamiento)
    for p in particles
        d = norm(p.r)
        delta = d + p.radius - radio_Recipiente

        if delta > 0.0
            n_wall = -p.r / d
            vn = dot(p.v, n_wall)

            Fn_mag = max(kn_wall * delta - gamma_wall * vn, 0.0)
            p.a += (Fn_mag * n_wall) / p.mass
        end
    end
end

# -- Funcion para aplicar las fuerzas de contacto entre las particulas -- #
function contacto_particulas!(particles::Vector{Particle{N, T}}) where {N, T}
    # Parametros de la interaccion entre particulas
    kn = 1.0e5
    gamma_n = 5.0e1
    num_p = length(particles)

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
                vn = dot(v_ij, n_ij)                # Componente normal de la velocidad relativa

                # Calculo de la fuerzas de contacto normales entre particulas
                # El modelo de fuerzas es linear spring-dashpot
                Fn_mag = max(kn * delta - gamma_n * vn, 0.0)
                Fn_vec = Fn_mag * n_ij
                
                # Aplicacion de la fuerza a las particulas i y j via 3era ley de Newton.
                particles[i].a += Fn_vec / p_i.mass
                particles[j].a -= Fn_vec / p_j.mass
            end
        end
    end
end

# -- Funcion que calcula todas las fuerzas que actuan sobre las particulas del sistema -- #
function fuerza_total!(particles::Vector{Particle{N, T}}, tiempo::T, radio_Recipiente::T) where {N, T}
    # Reiniciar aceleraciones
    for p in particles
        p.a = @SVector zeros(T, N)
    end

    # -- Aceleracion inercial asociado al movimiento de swirling del recipiente -- #
    excitacion_orbital_rampa!(particles, tiempo)

    # -- Fuerzas debidas al contacto con el recipiente circular -- #
    contenedor_circular!(particles, radio_Recipiente)

    # -- Fuerzas debidas al contacto entre particulas -- #
    contacto_particulas!(particles)
end

# -- Funcion para calcular la velocidad del marco de referencia del recipiente -- #
function velocidad_marco(tiempo::T) where {T}
    amplitud = 1.5
    frecuencia = 1.0
    tau = 1.0
    omega = 2.0 * pi * frecuencia

    factor_exp = exp(-tiempo/tau)
    
    A_t = amplitud * (1. - factor_exp * (1.0 + tiempo/tau))
    v_A = amplitud * (tiempo / tau^2) * factor_exp

    v_fx = v_A * cos(omega * tiempo) - A_t * omega * sin(omega * tiempo)
    v_fy = v_A * sin(omega * tiempo) + A_t * omega * cos(omega * tiempo)

    return @SVector [v_fx, v_fy]
end

# -- Funcion para calcular la energia cinetica total del sistema de particulas -- #
function calcular_energia_cinetica_absoluta!(particles::Vector{Particle{N, T}}, tiempo::T) where {N, T}
    V_f = velocidad_marco(tiempo)
    E_k = 0.0
    
    for p in particles
        v_abs = p.v + V_f
        E_k += 0.5 * p.mass * dot(v_abs, v_abs)
    end
    
    return E_k
end

# -- Funciones para exportar datos de la simulacion -- #

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

# -- Funcion para exportar datos extraidos de la simulacion -- #
function guardar_datos(archivo::String, tiempo::Float64, energia_cinetica::Float64, energia_potencial::Float64)
    open(archivo, "a") do io
        println(io, "$tiempo,$energia_cinetica,$energia_potencial,$(energia_cinetica + energia_potencial)")
    end
end

# -- Funcion principal para simular el sistema de particulas -- #
function simular_sistema()
    # -- Parametros del sistema -- #
    dimension_Sistema = 2
    n_Particulas = 1
    radio_Recipiente = 4.0          # En centimetros
    radio_Particula = 0.0           # En centimetros
    masa_Particula = 0.0            # En gramos

    # -- Inicializacion del sistema de particulas -- #
    sistema = [Particle(@SVector[0.0, 0.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5)]

    # -- Parametros de la ejecucion -- #
    dt = 1.0e-4
    tiempo_total = 1000.0 
    pasos = round(Int, tiempo_total / dt)
    frecuencia_guardado = 100
    
    # MALDITA SEA, ACUERDATE DE CAMBIAR EL NOMBRE DEL ARCHIVO DE SALIDA, NO TE VUELVAS A EQUIVOCAR
    archivo_salida_1 = "resultados/giro_swirling_150.xyz"
    archivo_salida_2 = "datos/datos_swirling_150.csv"
    open(archivo_salida_1, "w") do io end # Limpiar archivo si existe
    open(archivo_salida_2, "w") do io end # Limpiar archivo si existe

    println("Iniciando movimiento swirled...")
    
    for paso in 1:pasos
        velocity_verlet_step!(sistema, dt, fuerza_total!, dt * paso, radio_Recipiente)
        
        if paso % frecuencia_guardado == 0
            guardar_frame_xyz!(archivo_salida_1, sistema, paso * dt, radio_Recipiente)
            guardar_datos(archivo_salida_2, paso * dt, calcular_energia_cinetica_absoluta!(sistema, paso * dt), 0.0)
        end
    end
    
    println("¡Simulación terminada! Archivos generados: $archivo_salida_1 y $archivo_salida_2")
end

simular_sistema()