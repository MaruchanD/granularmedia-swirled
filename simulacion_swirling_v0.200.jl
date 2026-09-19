#=
Primer salto de versionado. A partir de aqui, las modificaciones seran de forma. El control de versionado
lo llevare con V.200-* donde el asterico refiere al numero del commit. Hacer mas commits.
=#
using StaticArrays
using LinearAlgebra

# == Estructura de datos para representar una particula == #
mutable struct Particle{N, T}
    # -- Propiedades de la particula -- #
    mass::T                 # Masa
    radius::T               # Radio
    inertia::T              # Momento de inercia

    # -- Grados de libertad traslacionales -- #
    r::SVector{N, T}        # Posición
    v::SVector{N, T}        # Velocidad
    a::SVector{N, T}        # Aceleración
    
    # -- Grados de libertad rotacionales (2D) -- #
    theta::T                # Ángulo de orientación
    omega::T                # Velocidad angular
    alpha::T                # Aceleración angular

    # -- Memoria cinemática (estado t) --#
    v_old::SVector{N, T}    # Velocidad en t
    omega_old::T            # Velocidad angular en t
    a_old::SVector{N, T}    # Aceleración en t
    alpha_old::T            # Aceleración angular en t

    # -- Memoria dinamica para balances de energía (pared y partículas) -- #
    f_pared::SVector{N, T}  # Fuerza contra la pared en t + dt
    tau_pared::T            # Torque contra la pared en t + dt
    f_pared_old::SVector{N, T}  # Fuerza contra la pared en t
    tau_pared_old::T        # Torque contra la pared en t
    
    f_part::SVector{N, T}   # Fuerza sobre la particula en t + dt
    tau_part::T             # Torque sobre la particula en t + dt
    f_part_old::SVector{N, T}   # Fuerza sobre la particula en t
    tau_part_old::T         # Torque sobre la particula en t
end

# == Funciones que generan las configuraciones de particulas. == #

# -- Funcion para generar el sistema de particulas iniciado a cero -- #
# Se ponen ciertos valores por defecto para la masa y el radio de las particulas, pero se pueden cambiar al llamar a la funcion.
function generar_sistema(N::Int, D::Int; m::Float64 = 1.0, R::Float64 = 1.0, theta::Float64 = 0.0, omega::Float64 = 0.0, alpha::Float64 = 0.0)
    # Se estan considerando por ahora que las particulas son discos uniformes, por lo que el momento de inercia es I = m*R^2/2.
    sistema = [Particle(
        #= Propiedades de la particula=#
        m,                          # Masa
        R,                          # Radio
        m*R^2/2,                    # Momento de Inercia (discos)
        #= Grados de libertad traslacional =#
        zeros(SVector{D, Float64}), # Posicion
        zeros(SVector{D, Float64}), # Velocidad
        zeros(SVector{D, Float64}), # Aceleracion
        #= Grados de libertad rotacional =#
        theta,                      # Angulo de orientacion
        omega,                      # Velocidad angular
        alpha,                      # Aceleracion angular
        #= Memoria cinematica =#
        zeros(SVector{D, Float64}), # Velocida en el paso t + dt
        0.0,                        # Velocidad angular en el paso t + dt
        zeros(SVector{D, Float64}), # Aceleración en el paso t
        0.0,                        # Aceleración angular en el paso t
        #= Memoria dinamica =#
        zeros(SVector{D, Float64}), # Fuerza contra la pared en t + dt
        0.0,                        # Torque contra la pared en t + dt
        zeros(SVector{D, Float64}), # Fuerza contra la pared en t
        0.0,                        # Torque contra la pared en t
        zeros(SVector{D, Float64}), # Fuerza debido a otra particula en t + dt
        0.0,                        # Torque debido a otra particula en t + dt
        zeros(SVector{D, Float64}), # Fuerza debido a otra particula en t
        0.0                         # Torque debido a otra particula en t
        ) for i in 1:N]
    println("El sistema de particulas de dimension ", D," con ", N, " particulas fue creado.")
    return sistema
end

# -- Funcion para generar una configuracion inicial de las particulas dentro de un recipiente circular -- #
# Se permite que las particulas se superpongan al inicio, pero luego se realiza un proceso de relajacion para minimizar las superposiciones.
# Ademas sea avisa al usuario si la fraccion de empaquetamiento supera el limite teorico de Random Close Packing en 2D (aprox 0.84).
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

# == Funciones que generan la dinamica del sistema de particulas == #

# -- Funcion para realizar un paso de integracion de Velocity Verlet -- #
# Se le pasa el arreglo que contiene las particulas, el paso de tiempo, la funcion que calcula las fuerzas, el tiempo actual,
# la tupla que contiene todos los parametros de la simulacion y el radio del recipiente. Se actualizan las posiciones, velocidades y aceleraciones de cada particula.
function velocity_verlet_step!(particles::Vector{Particle{N, T}}, dt::T, calc_forces!, tiempo::T, parametros_Generales::NamedTuple, energia_Mecanica::Vector{T}, radio_Recipiente::T) where {N, T}
    
    num_p = length(particles)

    # 1. Almacenar estados anteriores (estado t)
    for p in particles
        p.v_old = p.v
        p.omega_old = p.omega
        p.a_old = p.a
        p.alpha_old = p.alpha
        
        p.f_pared_old = p.f_pared
        p.tau_pared_old = p.tau_pared
        p.f_part_old = p.f_part
        p.tau_part_old = p.tau_part
    end
    # 2. Paso Predictor (estado t + dt estimado)
    for p in particles
        # -- Predicción de variables de traslacion -- #
        p.r = p.r + p.v * dt + 0.5 * p.a * dt^2
        p.v = p.v + p.a * dt  # Predicción con dt completo
        
        # -- Predicción de variables de rotacion -- #
        p.theta = mod(p.theta + p.omega * dt + 0.5 * p.alpha * dt^2, 2*pi)
        p.omega = p.omega + p.alpha * dt  # Predicción con dt completo
    end
    # 3. Evaluación de interacciones 
    # Se calculan las fuerzas usando las posiciones y velocidades predichas. 
    # Las variables p.a y p.alpha se sobrescriben con las nuevas aceleraciones reales (estado t + dt).
    calc_forces!(particles, tiempo + dt, parametros_Generales, radio_Recipiente)

    # 4. Paso Corrector (estado t + dt definitivo)
    for p in particles
        # -- Corrección de variables de traslacion -- #
        # Esto es v_new = v_pred + 0.5 * (a_new - a_old) * dt
        p.v = p.v + 0.5 * (p.a - p.a_old) * dt

        # -- Corrección de variables de rotacion -- #
        p.omega = p.omega + 0.5 * (p.alpha - p.alpha_old) * dt  # No sé explicar el signo menos
    end

    # 5. Calculo del Trabajo Inercial del paso
    trabajo_Inercial_step!(particles, dt, tiempo, parametros_Generales, energia_Mecanica)

    # 6. Calculo del Trabajo de la fuerza de contacto con la pared (y el torque)
    trabajo_Pared_step!(particles, dt, energia_Mecanica)

    # 7. Calculo del Trabajo de la fuerza de contacto entre particulas
    trabajo_Particulas_step!(particles, dt, energia_Mecanica)

    # 8. Calculo de la Energia Cinetica
    energia_Cinetica_step!(particles, tiempo, parametros_Generales, energia_Mecanica)

end

# -- Función pura para obtener el vector de aceleración inercial en cualquier instante -- #
function obtener_aceleracion_inercial(tiempo::T, (; amplitud, frecuencia, tau)) where {T}
    # Velocidad angular de la excitacion
    omega = 2.0 * pi * frecuencia

    # Envolvente para la amplitud y sus derivadas
    factor_exp = exp(-tiempo/tau)       
    
    A_t = amplitud * (1. - factor_exp * (1.0 + tiempo/tau))
    v_A = amplitud * (tiempo / tau^2) * factor_exp
    a_A = (amplitud / tau^2) * factor_exp * (1.0 - tiempo/tau)

    # Aceleracion inercial exacta
    aceleracion_x = -(a_A * cos(omega * tiempo) - 2.0 * v_A * omega * sin(omega * tiempo) - A_t * omega^2 * cos(omega * tiempo))
    aceleracion_y = -(a_A * sin(omega * tiempo) + 2.0 * v_A * omega * cos(omega * tiempo) - A_t * omega^2 * sin(omega * tiempo))
    
    return @SVector [aceleracion_x, aceleracion_y]
end

# -- Funcion para aplicar la excitacion orbital -- #
function excitacion_orbital_rampa!(particles::Vector{Particle{N,T}}, tiempo::T, parametros_Generales::NamedTuple) where {N, T}                         
    
    a_inercial = obtener_aceleracion_inercial(tiempo, parametros_Generales)
    
    # Aplicar la aceleracion inercial a cada particula
    for p in particles
        p.a += a_inercial
    end
end

# -- Funcion para aplicar las fuerzas de contacto entre las particulas y el contenedor circular -- #
function interaccion_recipiente!(particles::Vector{Particle{N, T}}, radio_Recipiente::T, (; dt, k_n_wall, gamma_n_wall, gamma_t_wall, mu_wall)) where {N, T}
    # El modelo de fuerzas usado es linear spring-dashpot (resorte + amortiguamiento)
    for p in particles
        d = norm(p.r)                               # Distancia de la particula al centro del recipiente
        delta = d + p.radius - radio_Recipiente     # Solapamiento de la particula con el recipiente

        # Aplicacion de las fuerzas de contacto normales y tangenciales entre particulas y contenedor
        if delta > 0.0
            n_wall = -p.r / d               # Direccion normal hacia el centro del recipiente
            r_c_wall = -p.radius * n_wall   # Vector desde el centro de la particula hasta el punto de contacto con la pared
            
            # Velocidad rotacional en el punto de contacto
            # En 2D, la velocidad rotacional es perpendicular al radio
            v_rot_p = p.omega * SVector(-r_c_wall[2], r_c_wall[1])

            # Velocidad total en el punto de contacto
            v = p.v + v_rot_p

            # Descomponemos la velocidad en el punto de contacto en las componentes normal y tangencial
            v_n_mag = dot(v, n_wall)        # Magnitud de la velocidad en la direccion normal a la pared del recipiente
            v_n_vec = v_n_mag * n_wall      # Componente normal de la velocidad

            v_t_vec = v - v_n_vec           # Componente tangencial de la velocidad
            v_t_mag = norm(v_t_vec)         # Magnitud de la velocidad en la direccion tangencial

            # Fuerza de contacto normal con la pared
            F_n_mag = max(k_n_wall * delta - gamma_n_wall * v_n_mag, 0.0)   # Esta comparacion evita que la fuerza sea negativa, es decir, evita que la fuerza sea atractiva.
            F_n = F_n_mag * n_wall

            F_contacto_pared = F_n

            # Fuerza de contacto tangencial con la pared
            if v_t_mag > 0.0
                # Vector Tangencial Unitario
                t_wall = v_t_vec / v_t_mag

                # Fuerza tangencial de prueba con un termino viscoso
                F_t_mag_prueba = gamma_t_wall * v_t_mag

                # Limite de Coulomb
                # La maxima fuerza permitida depende de la magnitud de la fuerza normal calculada
                F_t_mag_coulomb = mu_wall * F_n_mag

                # La magnitud final habra de ser el menor valor
                F_t_mag = min(F_t_mag_prueba, F_t_mag_coulomb)
                #= Se coloca el signo menos ya que la fuerza tangencial en este caso
                va siempre en el sentido contrario al vector t_wall ya que este se calcula en base
                a la componente tangencial de la velocidad en el punto de contacto de la particula
                con la pared =#
                F_t = -F_t_mag * t_wall
                F_contacto_pared += F_t

                #Calculo y aplicacion del torque debido a la fuerza tangencial
                tau = r_c_wall[1] * F_t[2] - r_c_wall[2] * F_t[1]
                p.alpha += tau / p.inertia
                p.tau_pared = tau
            end
            # Aplicar a la aceleración y registrar en la memoria de la partícula
            p.a += F_contacto_pared / p.mass
            p.f_pared = F_contacto_pared
        end
    end
end

# -- Funcion para aplicar las fuerzas de contacto entre las particulas -- #
function interaccion_particulas!(particles::Vector{Particle{N, T}}, (; dt, k_n, gamma_n, mu, gamma_t)) where {N, T}
    # Numero de particulas en el sistema
    num_p = length(particles)

    # Bucle doble para calcular las fuerzas de contacto entre todas las particulas del sistema
    for i in 1:num_p
        p_i = particles[i]
        for j in (i+1):num_p
            p_j = particles[j]

            # Definicion de la distancia entre particulas y la suma de sus radios.
            r_ij = p_i.r - p_j.r    # El vector que va de la particula j a la particula i
            d = norm(r_ij)          # Norma del vector de distancia entre las particulas
            suma_radios = p_i.radius + p_j.radius

            if d < suma_radios && d > 0.0
                # Vector Normal Unitario, que apunta de la particula j a la particula i
                n_ij = r_ij / d
                # Vector Tangencial Unitario
                t_ij = 

                # Vectores desde los centros de las particulas hasta el punto de contacto
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

                v_n_mag = norm(v_n)                 # Magnitud de la componente normal de la velocidad relativa
                v_t_mag = norm(v_t)                 # Magnitud de la componente tangencial de la velocidad relativa

                # Calculo de la fuerzas de contacto
                # Fuerza Normal
                # El modelo de fuerzas es linear spring-dashpot
                delta = suma_radios - d             # Solapamiento de las particulas
                F_n_mag = max(k_n * delta - gamma_n * v_n_mag, 0.0)
                F_n = F_n_mag * n_ij

                # Fuerza y torque total
                F_total_ij = F_n
                tau_i_total = 0.0
                tau_j_total = 0.0

                # Fuerza Tangencial
                if v_t_mag > 0.0
                    # Vector Tangencial Unitario
                    t_ij = v_t / v_t_mag
                    # Fuerza tangencial de prueba de tipo viscoso
                    F_t_mag_prueba = gamma_t * v_t_mag

                    # Limite de Coulomb
                    # La maxima fuerza permitida depende de la magnitud de la fuerza normal calculada
                    F_t_mag_coulomb = mu * F_n_mag

                    # La magnitud final es el menor de los dos valores
                    F_t_mag = min(F_t_mag_prueba, F_t_mag_coulomb)
                    F_t = -F_t_mag * t_ij                        # Revisar el signo de esta expresion.

                    # Calculo y aplicacion del torque debido a la fuerza tangencial
                    # Torque sobre la particula i
                    tau_i = r_ci[1] * F_t[2] - r_ci[2] * F_t[1]

                    # Torque sobre la particula j
                    # La fuerza tangencial sobre la particula j es opuesta a la de la particula i.
                    tau_j = r_cj[1]*(-F_t[2]) - r_cj[2]*(-F_t[1])

                    # Se suma la parte tangencial ahora
                    F_total_ij += F_t
                    tau_i_total = tau_i
                    tau_j_total = tau_j
                end
                # Aplicación en aceleración (Sobre i y j via 3era ley de Newton)
                particles[i].a += F_total_ij / p_i.mass
                particles[j].a -= F_total_ij / p_j.mass
                
                particles[i].alpha += tau_i_total / p_i.inertia
                particles[j].alpha += tau_j_total / p_j.inertia

                # Registro para el balance de energía
                particles[i].f_part += F_total_ij
                particles[j].f_part -= F_total_ij
                
                particles[i].tau_part += tau_i_total
                particles[j].tau_part += tau_j_total
            end
        end
    end
end

# -- Funcion que calcula todas las fuerzas que actuan sobre las particulas del sistema -- #
function fuerza_total!(particles::Vector{Particle{N, T}}, tiempo::T, parametros_Generales::NamedTuple, radio_Recipiente::T) where {N, T}
    # Reinicio de vectores
    for p in particles
        # -- Aceleraciones traslacionales -- #
        p.a = @SVector zeros(T, N)
        # -- Aceleraciones rotacionales -- #
        p.alpha = 0.0
        # -- Fuerzas del paso actual (t+dt) -- #
        p.f_pared = @SVector zeros(T, N)
        p.tau_pared = 0.0
        p.f_part = @SVector zeros(T, N)
        p.tau_part = 0.0
    end

    # -- Aceleracion inercial asociado al movimiento de swirling del recipiente -- #
    excitacion_orbital_rampa!(particles, tiempo, parametros_Generales)

    # -- Fuerzas debidas al contacto con el recipiente circular -- #
    interaccion_recipiente!(particles, radio_Recipiente, parametros_Generales)

    # -- Fuerzas debidas al contacto entre particulas -- #
    interaccion_particulas!(particles, parametros_Generales)
end

# == Funciones auxiliares para obtener informacion de la simulacion == #

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
function energia_Cinetica_step!(particles::Vector{Particle{N, T}}, tiempo::T, parametros_Generales::NamedTuple, energia_Mecanica::Vector{T}) where {N, T}
    # Primero se calcula la velocidad del marco de referencia del recipiente, y luego se calcula la energia cinetica total del sistema de particulas en el laboratorio.
    V_f = velocidad_marco(tiempo + parametros_Generales.dt, parametros_Generales)
    E_k_rel = 0.0
    E_k_traslacional = 0.0
    E_k_rotacional = 0.0
    
    for p in particles
        # Velocidad absoluta de la particual en el laboratorio
        v_abs = p.v + V_f

        # Energia cinetica relativa al marco giratorio
        E_k_rel += 0.5 * p.mass * dot(p.v, p.v)

        # Energia cinetica en el marco del laboratorio
        E_k_traslacional += 0.5 * p.mass * dot(v_abs, v_abs)
        E_k_rotacional += 0.5 * p.inertia * p.omega^2
    end
    # Energia cinetica traslacional y rotacional
    energia_Mecanica[1] = E_k_traslacional
    energia_Mecanica[2] = E_k_rotacional
    # Energia cinetica relativa al marco giratorio
    energia_Mecanica[3] = E_k_rel
    # Energia cinetica total
    energia_Mecanica[4] = E_k_traslacional + E_k_rotacional
end

# -- Funcion para calcular el trabajo inercial de la fuerza inercial despues de cada paso del integrador -- #
function trabajo_Inercial_step!(particles::Vector{Particle{N, T}}, dt::T, tiempo::T, parametros_Generales::NamedTuple, energia_Mecanica::Vector{T}) where {N, T}
    # Cálculo del Trabajo Inercial del paso
    a_inercial_old = obtener_aceleracion_inercial(tiempo, parametros_Generales)
    a_inercial_actual = obtener_aceleracion_inercial(tiempo + dt, parametros_Generales)
    a_inercial_media = 0.5 * (a_inercial_actual + a_inercial_old)
    potencia_inercial = 0.0
    
    for p in particles
        # Promedio exacto de la velocidad durante este intervalo dt
        v_media = 0.5 * (p.v + p.v_old) 
        potencia_inercial += dot(p.mass * a_inercial_media, v_media)
    end
    
    energia_Mecanica[7] += potencia_inercial * dt
end
# -- Funcion para calcular la energia disipada por choques con las paredes -- #
function trabajo_Pared_step!(particles::Vector{Particle{N, T}}, dt::T, energia_Mecanica::Vector{T}) where {N, T}
    # Variable para acumular trabajo dentro del bucle
    trabajo_pared_paso = 0.0

    for p in particles
        v_media = 0.5 * (p.v_old + p.v)
        omega_media = 0.5 * (p.omega_old + p.omega)
        f_pared_media = 0.5 * (p.f_pared_old + p.f_pared)
        tau_pared_media = 0.5 * (p.tau_pared_old + p.tau_pared)
        
        # Potencia mecánica de la pared (incluye tanto la elasticidad como la disipación)
        potencia_traslacional = dot(f_pared_media, v_media)
        potencia_rotacional = tau_pared_media * omega_media
        
        trabajo_pared_paso += (potencia_traslacional + potencia_rotacional) * dt

    end
    # Almacenar en el arreglo de energia
    energia_Mecanica[5] += trabajo_pared_paso
end

# -- Funcion para calcular la energia disipada por choques entre particulas -- #
function trabajo_Particulas_step!(particles::Vector{Particle{N, T}}, dt::T, energia_Mecanica::Vector{T}) where {N, T}
    # Variable para acumular el trabajo dentro del bucle
    trabajo_part_paso = 0.0
    for p in particles
        v_media = 0.5 * (p.v_old + p.v)
        omega_media = 0.5 * (p.omega_old + p.omega)
        f_part_media = 0.5 * (p.f_part_old + p.f_part)
        tau_part_media = 0.5 * (p.tau_part_old + p.tau_part)

        # Potencia mecanica entre particulas (elasticidad y disipacion)
        trabajo_part_paso += (dot(f_part_media, v_media) + tau_part_media * omega_media) * dt
    end
    # Almacenar en el arreglo de energía
    energia_Mecanica[6] += trabajo_part_paso
end

# -- Función para calcular la frecuencia de rotación del clúster -- #
function frecuencia_rotacion_cluster(particles::Vector{Particle{N, T}}) where {N, T}
    num_p = length(particles)
    
    # 1. Calcular el centro de masas y la velocidad del centro de masas
    r_cm = sum(p.r for p in particles) / num_p
    v_cm = sum(p.v for p in particles) / num_p

    sumatoria_f = 0.0
    for p in particles
        # Variables relativas al centro de masas
        r_star = p.r - r_cm
        v_star = p.v - v_cm
        
        r_star_sq = dot(r_star, r_star)
        
        if r_star_sq > 1e-12 # Evitar división por cero en partículas muy cerca del CM
            # Producto cruzado bidimensional: r_x * v_y - r_y * v_x
            cross_prod = r_star[1] * v_star[2] - r_star[2] * v_star[1]
            sumatoria_f += cross_prod / r_star_sq
        end
    end
    
    # La Ecuación 8 arroja la velocidad angular. Dividimos por 2*pi para obtener Hz.
    return (sumatoria_f / num_p) / (2 * pi)
end

# -- Función para calcular el momento de inercia respecto al eje del recipiente -- #
function inercia_cluster(particles::Vector{Particle{N, T}}) where {N, T}
    I_total = 0.0
    for p in particles
        # Distancia al cuadrado respecto al centro del contenedor (origen)
        r_sq = dot(p.r, p.r) 
        # Teorema de ejes paralelos: I_eje = I_cm + m * d^2
        I_total += p.inertia + p.mass * r_sq
    end
    return I_total
end

# == Funciones para exportar datos de la simulacion == #

# Exportacion a OVITO
function guardar_frame_xyz(archivo::String, particles::Vector{Particle{N, T}}, tiempo::Float64, radio_Recipiente::T) where {N, T}
    open(archivo, "a") do io
        # -- Número de partículas + contenedor -- #
        println(io, length(particles)+1)

        # -- Propiedades del archivo para ser leidas por el programa OVITO -- #
        println(io, "Properties=species:S:1:pos:R:3:radius:R:1:Theta:R:1:OmegaZ:R:1 Time=$tiempo")

        # -- Variables de las particulas -- #
        for p in particles
            # Se agregan las propiedades de rotacion (Theta y OmegaZ) al archivo de salida
            println(io, "Granulo $(p.r[1]) $(p.r[2]) 0.0 $(p.radius) $(p.theta) $(p.omega)")
        end

        # -- Variables del recipiente -- #
        println(io, "Contenedor 0.0 0.0 0.0 $radio_Recipiente 0.0 0.0")
    end
end

function guardar_datos(archivo::String, tiempo::Float64, energia_Mecanica::Vector{T}) where {T}
    energia_Cinetica_r = energia_Mecanica[2]              
    energia_Cinetica_relativa = energia_Mecanica[3]     
    trabajo_Pared = energia_Mecanica[5]
    trabajo_Particulas = energia_Mecanica[6]
    trabajo_Inercial = energia_Mecanica[7]

    # K_total = K_traslacional_relativa + K_rotacional
    energia_Cinetica_Total = energia_Cinetica_relativa + energia_Cinetica_r
    
    # Teorema del Trabajo y la Energía: ΔK = W_neto
    # Como la partícula parte del reposo, K_inicial = 0.0
    balance_diagnostico = energia_Cinetica_Total - trabajo_Inercial - trabajo_Pared - trabajo_Particulas

    open(archivo, "a") do io
        # Imprime los datos que consideres, pero observa puntualmente balance_diagnostico
        println(io, "$tiempo,$energia_Cinetica_relativa,$energia_Cinetica_r,$trabajo_Inercial,$trabajo_Pared,$trabajo_Particulas,$balance_diagnostico")
    end
end

# -- Función para guardar las variables macroscópicas del clúster -- #
function guardar_datos_macroscopicos(archivo::String, tiempo::Float64, num_p::Int, f_rot::Float64, I_cluster::Float64)
    open(archivo, "a") do io
        println(io, "$tiempo,$num_p,$f_rot,$I_cluster")
    end
end

# -- Funcion principal para simular el sistema de particulas -- #
function simular_sistema(N_input::Union{Int, Nothing} = nothing)
    # Condicional para verificar si se paso un numero en la llamada de la funcion
    numero_Particulas = if isnothing(N_input)
        # -- Solicitud dinámica del número de partículas al usuario -- #
        print("Introduzca el número de partículas del sistema: ")
        numero_Particulas = parse(Int, readline())
    else
        N_input
    end

    # -- Parametros del sistema -- #
    dimension_Sistema = 2
    radio_Recipiente = 6.0          # En centimetros
    radio_Particula = 0.5           # En centimetros
    masa_Particula = 1.0            # En gramos

    # -- Parametros para la dinamica del sistema -- #
    # LAS COMAS SON IMPORTANTES PORQUE SON ELEMENTOS DE UNA TUPLA, SI NO SE PONEN LAS COMAS, EL COMPILADOR NO LOS RECONOCE COMO ELEMENTOS DE LA TUPLA.
    parametros_Generales = (
        # -- Parametro del paso de simulacion -- #
        dt = 1.0e-4,

        # -- Parametros de la excitacion orbital -- #
        amplitud = 1.5,     # Radio de la excitacion (1.0 a 5.0 cm)
        frecuencia = 1.0,   # Frecuencia de la excitacion (0.1 - 5.0 Hz)
        tau = 0.5,          # Tiempo de rampa (s)

        # -- Parametros de colisiones particula-recipiente -- #
        
        # -- Normales -- #
        k_n_wall = 7.0e5, #
        gamma_n_wall = 1.0e3, #

        # -- Tangenciales -- #
        gamma_t_wall = 1.0e3, #
        mu_wall = 0.5, #

        # -- Parametros de colisiones particula-particula -- #
        
        # -- Normales -- #
        k_n = 7.0e5, #
        gamma_n = 1.0e3, #

        # -- Tangenciales -- #
        gamma_t = 1.0e3, #
        mu = 0.5 #
    )

    # -- Inicializacion del sistema de particulas -- #
    
    sistema = generar_sistema(numero_Particulas, dimension_Sistema, R = radio_Particula, m = masa_Particula)
    generar_configuracion!(sistema, radio_Recipiente)
    
    # Sistema de particulas de prueba para verificar funciones o estabilidad de la simulacion
    #=
    sistema = [
        Particle(@SVector[5.0, 0.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5, 0.0, 0.0, 0.0, 0.5*masa_Particula*radio_Particula^2, zeros(SVector{2, Float64}), 0.0)#=, 
        Particle(@SVector[-5.0, 0.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5, 0.0, 0.0, 0.0, 0.5*masa_Particula*radio_Particula^2, zeros(SVector{2, Float64}), 0.0)=#
    ]
    =#
    # -- Parametros de la ejecucion -- #
    dt = parametros_Generales.dt                    # En segundos (s)
    tiempo_total = 50.0                             # En segundos (s)
    pasos = round(Int, tiempo_total / dt)
    frecuencia_guardado = 100                       # Cada 0.01 segundos (s) se guardan los datos de la simulacion

    # -- Archivos de salida -- #
    amp = parametros_Generales.amplitud
    frec = parametros_Generales.frecuencia
    # Nombre de los archivos de salida para los datos de la simulacion y para los datos extraidos de la simulacion
    # MALDITA SEA, ACUERDATE DE CAMBIAR EL NOMBRE DEL ARCHIVO DE SALIDA, NO TE VUELVAS A EQUIVOCAR
    archivo_salida_1 = "resultados/simulacion_N$(numero_Particulas)_A$(amp)_F$(frec).xyz"
    archivo_salida_2 = "datos/datos_energia_N$(numero_Particulas)_A$(amp)_F$(frec).csv"
    archivo_salida_3 = "datos/datos_macroscopicos_N$(numero_Particulas)_A$(amp)_F$(frec).csv"

    # Apertura y limpieza de los archivos de salida
    open(archivo_salida_1, "w") do io end
    open(archivo_salida_2, "w") do io end
    open(archivo_salida_2, "w") do io end

    # -- Arreglo para registrar la energia del sistema -- # 
    energia_Mecanica = [
        0.0,    # [1] Energia cinetica traslacional real
        0.0,    # [2] Energia cinetica rotacional
        0.0,    # [3] Energia cinetica relativa
        0.0,    # [4] Energia cinetica total real

        # Variables energia debidas al trabajo de las fuerzas de contacto
        0.0,    # [5] Trabajo Pared acumulado
        0.0,    # [6] Trabajo particulas acumulado

        # Variable de energia debido al trabajo de la fuerza inercial
        0.0     # [7] Trabajo inercial
        ]

    println("Iniciando movimiento swirled...")
    
    # -- Bucle principal de la simulacion -- #
    for paso in 0:pasos
        # -- Integracion de las ecuaciones de movimiento mediante el metodo de Velocity Verlet -- #
        velocity_verlet_step!(sistema, dt, fuerza_total!, dt * paso, parametros_Generales, energia_Mecanica, radio_Recipiente)
        
        # -- Guardado de los datos de la simulacion y de los datos extraidos de la simulacion -- #
        if paso % frecuencia_guardado == 0
            tiempo_actual = paso * dt

            # 1. Cálculos macroscópicos del clúster
            f_rot = frecuencia_rotacion_cluster(sistema)
            I_cluster = inercia_cluster(sistema)

            # 2. Guardado en archivos
            guardar_frame_xyz(archivo_salida_1, sistema, tiempo_actual, radio_Recipiente)
            guardar_datos(archivo_salida_2, tiempo_actual, energia_Mecanica)
            guardar_datos_macroscopicos(archivo_salida_3, tiempo_actual, numero_Particulas, f_rot, I_cluster)
        end
    end
    
    println("¡Simulación terminada!\nArchivos generados:\n- $archivo_salida_1\n- $archivo_salida_2\n- $archivo_salida_3")
end

#simular_sistema()