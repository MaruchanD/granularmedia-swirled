# Las versiones 0.14* van a considerar el movimiento de rotacion alrededor del eje de simetria de las particulas.
# Se hacen modificaciones a la estructura de las particulas, al integrador y a las funciones de fuerzas para considerar el torque. 
# En esta primera version, solo se modifica la estructura y el integrador, dejando las funciones de fuerza sin modificar. 
# Esto con el objetivo de probar la implementacion del torque y la rotacion de las particulas.
# Tambien se modifica la funcion que guarda los frames, se obvia el generador de configuraciones. Se probaran configuraciones ya generadas.
# Se agrega una ligera modificacion al integrador para poder pasar el tamaño del recipiente.

using StaticArrays
using LinearAlgebra

mutable struct Particle{N, T}
    # -- Grados de liberta traslacionales -- #
    r::SVector{N, T}      # Posición
    v::SVector{N, T}      # Velocidad
    a::SVector{N, T}      # Aceleración
    mass::T               # Masa
    radius::T             # Radio

    # -- Grados de libertad rotacionales (2D) -- #
    theta::T              # Ángulo de orientación (radianes)
    omega::T              # Velocidad angular (radianes/s)
    alpha::T              # Aceleración angular (radianes/s^2)
    inertia::T            # Momento de inercia 
end

function velocity_verlet_step!(particles::Vector{Particle{N, T}}, dt::T, calc_forces!, tiempo::T, radio_Recipiente::T) where {N, T}
    for p in particles
        # -- Calculo de variables de traslacion --#
        p.r = p.r + p.v * dt + 0.5 * p.a * dt^2
        p.v = p.v + 0.5 * p.a * dt

        # -- Calculo de variables de rotacion --#
        p.theta = mod(p.theta + p.omega * dt + 0.5 * p.alpha * dt^2, 2*pi)  # Mantener el ángulo en el rango [0, 2π)
        p.omega = p.omega + 0.5 * p.alpha * dt
    end
    
    calc_forces!(particles, tiempo + dt, radio_Recipiente)
    
    for p in particles
        #-- Calculo de variables de traslacion --#
        p.v = p.v + 0.5 * p.a * dt

        #-- Calculo de variables de rotacion --#
        p.omega = p.omega + 0.5 * p.alpha * dt
    end
end

function excitacion_orbital_rampa!(particles::Vector{Particle{N,T}}, tiempo::T) where {N, T}
    # Parametros para la excitacion
    amplitud = 3.0                     # Radio de la excitacion (1.0 a 5.0 cm)
    frecuencia = 2.5                    # Frecuencia de la excitacion (0.1 - 5.0 Hz)
    tau = 1.0                           # Tiempo de rampa (s)
    omega = 2.0 * pi * frecuencia       # Facilidad para construir las cuentas

    # Envolvente para la amplitud y sus derivadas
    factor_exp = exp(-tiempo/tau)       #Factor exponencial
    
    A_t = amplitud * (1. - factor_exp * (1.0 + tiempo/tau))
    v_A = amplitud * (tiempo / tau^2) * factor_exp
    a_A = (amplitud / tau^2) * factor_exp * (1.0 - tiempo/tau)

    # Aceleracion inercial exacta
    aceleracion_x = (a_A * cos(omega * tiempo) - 2.0 * v_A * (omega) * sin(omega * tiempo) - A_t * (omega)^2 * cos(omega * tiempo))
    aceleracion_y = (a_A * sin(omega * tiempo) + 2.0 * v_A * (omega) * cos(omega * tiempo) - A_t * (omega)^2 * sin(omega * tiempo))
    aceleracion_inercial = @SVector [aceleracion_x, aceleracion_y]

    # Aplicar la aceleracion inercial a cada particula
    for p in particles
        p.a += aceleracion_inercial
    end
end

function contenedor_circular!(particles::Vector{Particle{N, T}}, radio_Recipiente::T) where {N, T}
    # -- Parametros de la interaccion con el contenedor --#
    # -- Normales --#
    kn_wall = 1.0e5
    gamma_wall = 5.0e3

    # -- Tangenciales --#
    mu_pared = 0.3
    gamma_t_pared = 8.0e2

    for p in particles
        d = norm(p.r)
        delta = d + p.radius - radio_Recipiente

        # Aplicacion de las fuerzas de contacto
        if delta > 0.0
            n_wall = -p.r / d               # Direccion normal al centro
            r_c_wall = -p.radius * n_wall   # Vector desde el centro de la particula hasta el punto de contacto con la pared

            # Velocidad rotacional en el punto de contacto
            v_rot_p = p.omega * SVector(-r_c_wall[2], r_c_wall[1])  # En 2D, la velocidad rotacional es perpendicular al radio

            v = p.v + v_rot_p # Velocidad total en el punto de contacto
            
            # Descomposicion de la velocidad en el punto de contacto en componentes normal y tangencial
            v_n_mag = dot(v, n_wall)        # Magnitud de la velocidad en la direccion normal
            v_n_vec = v_n_mag * n_wall      # Componente normal de la velocidad
            
            v_t_vec = v - v_n_vec           # Componente tangencial de la velocidad
            v_t_mag = norm(v_t_vec)         # Magnitud de la velocidad en la direccion tangencial

            # Fuerza de contacto normal con la pared
            F_n_mag = max(kn_wall * delta - gamma_wall * v_n_mag, 0.0)
            F_n_vec = F_n_mag * n_wall

            # Aplicacion de la fuerza normal
            p.a += F_n_vec / p.mass
            
            # Fuerza de contacto tangencial con la pared
            if v_t_mag > 0.0
                # Vector Tangencial Unitario
                t_wall = v_t_vec / v_t_mag

                # Fuerza tangencial de prueba (viscosa)
                Ft_prueba = gamma_t_pared * v_t_mag

                # Limite de Coulomb
                # La maxima fuerza permitida depende de la magnitud de la fuerza normal calculada
                Ft_limite_coulomb = mu_pared * F_n_mag

                # La magnitud final es el menor de los dos valores
                F_t_mag = min(Ft_prueba, Ft_limite_coulomb)
                F_t_vec = -F_t_mag * t_wall

                # Aplicacion de la fuerza tangencial lineal
                p.a += F_t_vec / p.mass

                # Calculo y aplicacion del torque debido a la fuerza tangencial
                tau = r_c_wall[1]*F_t_vec[2] - r_c_wall[2]*F_t_vec[1]
                p.alpha += tau / p.inertia
            end
        end
    end
end

function contacto_particulas!(particles::Vector{Particle{N, T}}) where {N, T}
    # -- Parametros de la interaccion entre particulas -- #
    # -- Normales -- #
    k_n = 1.0e5
    gamma_n = 5.0e1

    # -- Tangenciales -- #
    mu = 0.5
    gamma_t = 1.0e3

    num_p = length(particles)

    for i in 1:num_p
        p_i = particles[i]
        for j in (1+i):num_p
            p_j = particles[j]

            r_ij = p_i.r - p_j.r    # El vector que va de la particula j a la particula i
            d = norm(r_ij)          # Norma del vector de distancia entre las particulas
            suma_radios = p_i.radius + p_j.radius

            if d < suma_radios && d > 0.0
                # Vector Normal Unitario, apunta de la particula j a la particula i
                n_ij = r_ij / d

                # Vectores desde los centros hasta el punto de contacto
                r_ci = -p_i.radius * n_ij  # Vector desde el centro de la particula i hasta el punto de contacto
                r_cj = p_j.radius * n_ij   # Vector desde el centro de la particula j hasta el punto de contacto

                # Velocidades tangenciales debido a la rotacion (producto cruzado en 2D)
                v_rot_i = p_i.omega * SVector(-r_ci[2], r_ci[1])  # Velocidad rotacional de la particula i en el punto de contacto
                v_rot_j = p_j.omega * SVector(-r_cj[2], r_cj[1])  # Velocidad rotacional de la particula j en el punto de contacto

                # Velocidades totales en el punto de contacto
                v_ci = p_i.v + v_rot_i
                v_cj = p_j.v + v_rot_j

                # Velocidad relativa en el punto de contacto
                v_ij = v_ci - v_cj

                # Descomposicion de la velocidad relativa en componentes normal y tangencial
                v_n = dot(v_ij, n_ij) * n_ij        # Componenten normal de la velocidad relativa
                v_t = v_ij - v_n                    # Componente tangencial de la velocidad relativa

                v_n_mag = dot(v_ij, n_ij)           # Magnitud de la componente normal de la velocidad relativa
                v_t_mag = norm(v_t)                 # Magnitud de la componente tangencial de la velocidad relativa

                # Calculo de las fuerzas de contacto
                # Fuerza Normal
                delta = suma_radios - d                     # Solapamiento de las particulas
                Fn_mag = k_n * delta - gamma_n * v_n_mag
                Fn_vec = max(Fn_mag, 0.0) * n_ij

                # Aplicacion de la fuerza normal
                particles[i].a += Fn_vec / p_i.mass
                particles[j].a -= Fn_vec / p_j.mass

                # Fuerza Tangencial
                if v_t_mag > 0.0
                    # Vector Tangencial Unitario
                    t_ij = v_t / v_t_mag

                    # Fuerza tangencial de prueba (viscosa)
                    Ft_prueba = gamma_t * v_t_mag

                    # Limite de Coulomb
                    # La maxima fuerza permitida depende de la magnitud de la fuerza normal calculada
                    F_limite_coulomb = mu * Fn_mag

                    # La magnitud final es el menor de los dos valores
                    Ft_mag = min(Ft_prueba, F_limite_coulomb)
                    Ft_vec = -Ft_mag * t_ij

                    # Aplicacion de la fuerza tangencial
                    particles[i].a += Ft_vec / p_i.mass
                    particles[j].a -= Ft_vec / p_j.mass

                    # Calculo y aplicacion del torque debido a la fuerza tangencial
                    # Torque sobre la particula i
                    tau_i = r_ci[1]*Ft_vec[2] - r_ci[2]*Ft_vec[1]
                    particles[i].alpha += tau_i / p_i.inertia

                    # Torque sobre la particula j
                    # La fuerza tangencial sobre la particula j es opuesta a la de la particula i.
                    tau_j = r_cj[1]*(-Ft_vec[2]) - r_cj[2]*(-Ft_vec[1])
                    particles[j].alpha += tau_j / p_j.inertia
                end
            end
        end
    end
end

function fuerza_total!(particles::Vector{Particle{N, T}}, tiempo::T, radio_Recipiente::T) where {N, T}
    # Reiniciar aceleraciones
    for p in particles
        # -- Aceleraciones traslacionales -- #
        p.a = @SVector zeros(T, N)
        # -- Aceleraciones rotacionales -- #
        p.alpha = 0.0
    end

    # Aplicar todas las fuerzas
    excitacion_orbital_rampa!(particles, tiempo)
    contenedor_circular!(particles, radio_Recipiente)
    contacto_particulas!(particles)
end

# Exportacion a OVITO
function guardar_frame_xyz!(archivo::String, particles::Vector{Particle{N, T}}, tiempo::Float64, radio_Recipiente::T) where {N, T}
    open(archivo, "a") do io
        println(io, length(particles)+1)  # Número de partículas + contenedor
        println(io, "Properties=species:S:1:pos:R:3:radius:R:1:Theta:R:1:OmegaZ:R:1 Time=$tiempo")
        for p in particles
            # Se agregan las propiedades de rotacion (Theta y OmegaZ) al archivo de salida
            println(io, "Granulo $(p.r[1]) $(p.r[2]) 0.0 $(p.radius) $(p.theta) $(p.omega)")
        end
        println(io, "Contenedor 0.0 0.0 0.0 $radio_Recipiente 0.0 0.0")
    end
end

function simular_sistema()
    sistema = [
        Particle(@SVector[0.0, 0.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5, 0.0, 1.0, 0.0, 0.5),
        Particle(@SVector[2.0, 1.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5, 2*pi, 0.0, 0.0, 0.5),
        Particle(@SVector[-2.5, -1.5], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5, 0.0, 3.0, 0.0, 0.5)
    ]
    # Parametros
    radio_Recipiente = 5.0

    dt = 1.0e-4
    tiempo_total = 10.0 
    pasos = round(Int, tiempo_total / dt)
    frecuencia_guardado = 100 
    
    archivo_salida = "resultados/giro_swirling_140.xyz" # MALDITA SEA, ACUERDATE DE CAMBIAR EL NOMBRE DEL ARCHIVO DE SALIDA, NO TE VUELVAS A EQUIVOCAR
    open(archivo_salida, "w") do io end # Limpiar archivo si existe

    println("Iniciando movimiento swirled...")
    
    for paso in 1:pasos
        velocity_verlet_step!(sistema, dt, fuerza_total!, dt * paso, radio_Recipiente)
        
        if paso % frecuencia_guardado == 0
            guardar_frame_xyz!(archivo_salida, sistema, paso * dt, radio_Recipiente)
        end
    end
    
    println("¡Simulación terminada! Archivo generado: $archivo_salida")
end

simular_sistema()