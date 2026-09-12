using StaticArrays
using LinearAlgebra

# == Estructura de datos para representar una particula == #
mutable struct Particle{N, T}
    # -- Grados de liberta traslacionales -- #
    r::SVector{N, T}      # Posición
    v::SVector{N, T}      # Velocidad
    a::SVector{N, T}      # Aceleración
    mass::T               # Masa
    radius::T             # Radio

    # -- Grados de libertad rotacionales (2D) -- #
    theta::T              # Ángulo de orientación
    omega::T              # Velocidad angular
    alpha::T              # Aceleración angular
    inertia::T            # Momento de inercia 
end

# == Funciones que generan las configuraciones de particulas. == #

# -- Funcion para generar el sistema de particulas iniciado a cero -- #
# Se ponen ciertos valores por defecto para la masa y el radio de las particulas, pero se pueden cambiar al llamar a la funcion.
function generar_sistema(N::Int, D::Int; m::Float64 = 1.0, R::Float64 = 1.0, theta::Float64 = 0.0, omega::Float64 = 0.0, alpha::Float64 = 0.0)
    # Se estan considerando por ahora que las particulas son discos uniformes, por lo que el momento de inercia es I = m*R^2/2.
    sistema = [Particle(zeros(SVector{D, Float64}), zeros(SVector{D, Float64}), zeros(SVector{D, Float64}), m, R, theta, omega, alpha, m*R^2/2) for i in 1:N]
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
    a_old = [p.a for p in particles] # Se guarda la aceleracion del paso anterior
    alpha_old = [p.alpha for p in particles] # Se guarda la aceleracion angular del paso anterior
    v_old = [p.v for p in particles] # Se guarda la velocidad real al inicio del paso
    omega_old = [p.omega for p in particles] # Se guarda la velocidad angular al inicio del paso

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
    calc_forces!(particles, tiempo + dt, parametros_Generales, energia_Mecanica, radio_Recipiente)

    # 4. Paso Corrector (estado t + dt definitivo)
    for i in 1:num_p
        p = particles[i]
        
        # -- Corrección de variables de traslacion -- #
        # Esto es v_new = v_pred + 0.5 * (a_new - a_old) * dt
        p.v = p.v + 0.5 * (p.a - a_old[i]) * dt

        # -- Corrección de variables de rotacion -- #
        p.omega = p.omega + 0.5 * (p.alpha - alpha_old[i]) * dt
    end

    # 5. Calculo del Trabajo Inercial del paso
    trabajo_Inercial_step!(particles, dt, tiempo, parametros_Generales, v_old, energia_Mecanica)

    # 6. Calculo de la Potencia y Energia Disipada pared paso
    energia_Disipada_pared_step!(particles, v_old, omega_old, parametros_Generales, energia_Mecanica, radio_Recipiente)

    # 7. Calculo de la Energia Cinetica
    energia_Cinetica_step!(particles, tiempo, v_old, omega_old, parametros_Generales, energia_Mecanica)

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
function contenedor_circular!(particles::Vector{Particle{N, T}}, radio_Recipiente::T, (; k_n_wall, gamma_n_wall, gamma_t_wall, mu_wall), energia_Mecanica::Vector{T}) where {N, T}
    # Variable para acumular la energia potencial particula-pared
    energia_Potencial_part_pared = 0.0

    # Arreglo que contiene las fuerzas de contacto contra el contenedor de cada una de las particulas
    F_c = SVector{N, T}[]

    # El modelo de fuerzas usado es linear spring-dashpot (resorte + amortiguamiento)
    for p in particles
        d = norm(p.r)                               # Distancia de la particula al centro del recipiente
        delta = d + p.radius - radio_Recipiente     # Solapamiento de la particula con el recipiente

        # Aplicacion de las fuerzas de contacto normales y tangenciales entre particulas y contenedor
        if delta > 0.0
            n_wall = -p.r / d                       # Direccion normal hacia el centro del recipiente
            t_wall = SVector(-n_wall[2], n_wall[1]) # Direccion tangencial
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

            if v_t_mag > 0.0
                # Fuerza de contacto tangencial con la pared
                # Fuerza tangencial de prueba con un termino viscoso
                F_t_mag_prueba = gamma_t_wall * v_t_mag
                # Limite de Coulomb
                # La minima fuerza permitida depende de la magnitud de la fuerza normal calculada
                F_t_mag_coulomb = mu_wall * F_n_mag

                # La magnitud final habra de ser el menor valor
                F_t_mag = min(F_t_mag_prueba, F_t_mag_coulomb)
                #= Se coloca el signo menos ya que la fuerza tangencial en este caso
                va siempre en el sentido contrario al vector t_wall ya que este se calcula en base
                a la componente tangencial de la velocidad en el punto de contacto de la particula
                con la pared =#
                F_t = -F_t_mag * t_wall

                #Calculo y aplicacion del torque debido a la fuerza tangencial
                tau = r_c_wall[1] * F_t[2] - r_c_wall[2] * F_t[1]
                p.alpha += tau / p.inertia
            end
            push!(F_c, F_n + F_t)
            # Calculo de la energia potencial de la interaccion particula-pared
            energia_Potencial_part_pared += 0.5 * k_n_wall * delta^2
        end
    end
    # Registro de la energia potencial particula-pared un instante de tiempo.
    energia_Mecanica[4] = energia_Potencial_part_pared

    return F_c
end

# -- Funcion para aplicar las fuerzas de contacto entre las particulas -- #
function contacto_particulas!(particles::Vector{Particle{N, T}}, (; dt, k_n, gamma_n, mu, gamma_t), energia_Mecanica::Vector{T}) where {N, T}
    # Variable para acumular la energia potencial particula-particula
    energia_Potencial_part_part = 0.0
    # Variable para acumular la potencia disipada particula-particula
    potencia_Disipada_part_part = 0.0

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

                # Aplicacion de la fuerza normal (por tercera ley de Newton)
                particles[i].a += F_n / p_i.mass
                particles[j].a -= F_n / p_j.mass

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
                    F_t_vec = -F_t_mag * t_ij                        # Revisar el signo de esta expresion.

                    # Aplicacion de la fuerza tangencial
                    particles[i].a += F_t_vec / p_i.mass
                    particles[j].a -= F_t_vec / p_j.mass

                    # Calculo y aplicacion del torque debido a la fuerza tangencial
                    # Torque sobre la particula i
                    tau_i = r_ci[1] * F_t_vec[2] - r_ci[2] * F_t_vec[1]
                    particles[i].alpha += tau_i / p_i.inertia

                    # Torque sobre la particula j
                    # La fuerza tangencial sobre la particula j es opuesta a la de la particula i.
                    tau_j = r_cj[1]*(-F_t_vec[2]) - r_cj[2]*(-F_t_vec[1])
                    particles[j].alpha += tau_j / p_j.inertia
  
                    # Calculo de la potencia disipada en la interaccion particula-particula
                    # Aqui se escribe la potencia disipada (F*v) en 2 partes, las cuales corresponden
                    # a los dos terminos viscosos que se estan considerando: el normal con gamma_n
                    # y el tangencial que viene representado por el resultado de F_t_mag.
                    potencia_Disipada_part_part += gamma_n * v_n_mag^2 + F_t_mag * v_t_mag

                    # ESTE IF NECESITA ALGUNAS PRUEBAS ADICIONALES PARA SABER SI ES CONFIABLE
                end

                # Calculo de la energia potencial de la interaccion particula-particula
                energia_Potencial_part_part += 0.5 * k_n * delta^2
            end
        end
    end
    # Registro de la energia potencial y disipada particula-particula en un instante de tiempo.
    energia_Mecanica[6] = energia_Potencial_part_part
    energia_Mecanica[7] += potencia_Disipada_part_part*dt
end

# -- Funcion que calcula todas las fuerzas que actuan sobre las particulas del sistema -- #
function fuerza_total!(particles::Vector{Particle{N, T}}, tiempo::T, parametros_Generales::NamedTuple, energia_Mecanica::Vector{T}, radio_Recipiente::T) where {N, T}
    # Reiniciar aceleraciones
    for p in particles
        p.a = @SVector zeros(T, N)
        # -- Aceleraciones rotacionales -- #
        p.alpha = 0.0
    end

    # -- Aceleracion inercial asociado al movimiento de swirling del recipiente -- #
    excitacion_orbital_rampa!(particles, tiempo, parametros_Generales)

    # -- Fuerzas debidas al contacto con el recipiente circular -- #
    contenedor_circular!(particles, radio_Recipiente, parametros_Generales, energia_Mecanica)

    # -- Fuerzas debidas al contacto entre particulas -- #
    contacto_particulas!(particles, parametros_Generales, energia_Mecanica)
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
    
    for (i,p) in enumerate(particles)
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

    # Energia cinetica total
    energia_Mecanica[3] = E_k_traslacional + E_k_rotacional
    
    # Energia cinetica relativa al marco giratorio
    energia_Mecanica[9] = E_k_rel
end

# -- Funcion para calcular el trabajo inercial de la fuerza inercial despues de cada paso del integrador -- #
function trabajo_Inercial_step!(particles::Vector{Particle{N, T}}, dt::T, tiempo::T, parametros_Generales::NamedTuple, v_old::Vector{SVector{N, T}}, energia_Mecanica::Vector{T}) where {N, T}
    # Cálculo del Trabajo Inercial del paso
    a_inercial_old = obtener_aceleracion_inercial(tiempo, parametros_Generales)
    a_inercial_actual = obtener_aceleracion_inercial(tiempo + dt, parametros_Generales)
    a_inercial_media = 0.5 * (a_inercial_actual + a_inercial_old)
    potencia_inercial = 0.0
    
    for (i, p) in enumerate(particles)
        # Promedio exacto de la velocidad durante este intervalo dt
        v_media = 0.5 * (v_old[i] + p.v) 
        potencia_inercial += dot(p.mass * a_inercial_media, v_media)
    end
    
    energia_Mecanica[8] += potencia_inercial * dt
end
# -- Funcion para calcular la energia disipada por choques con las paredes -- #
function energia_Disipada_pared_step!(particles::Vector{Particle{N, T}}, v_old::Vector{SVector{N, T}}, omega_old::Vector{T}, (; dt, k_n_wall, gamma_n_wall, gamma_t_wall, mu_wall), energia_Mecanica::Vector{T}, radio_Recipiente::T) where {N, T}
    potencia_disipada_pared = 0.0

    for (i, p) in enumerate(particles)
        d = norm(p.r)
        delta = d + p.radius - radio_Recipiente

        if delta > 0.0
            n_wall = -p.r / d
            r_c_wall = -p.radius * n_wall

            #=
            # Velocidad media en el intervalo [t, t + dt]
            v_tras_media = 0.5 * (v_old[i] + p.v)
            omega_media = 0.5 * (omega_old[i] + p.omega)
            v_rot_media = omega_media * SVector(-r_c_wall[2], r_c_wall[1])
            v_cont_media = v_tras_media + v_rot_media

            # Componentes de la velocidad media
            v_n_media_mag = dot(v_cont_media, n_wall)
            v_t_media_vec = v_cont_media - v_n_media_mag * n_wall
            v_t_media_mag = norm(v_t_media_vec)=#

            # Calculo de las componentes de la velocidad en el paso t
            v_rot_old = omega_old[i] * SVector(-r_c_wall[2], r_c_wall[1])
            v_cont_old = v_old[i] + v_rot_old
            v_n_old_mag = dot(v_cont_old, n_wall)
            v_t_old = v_cont_old - v_n_old_mag * n_wall
            v_t_old_mag = norm(v_t_old)

            # Calculo de las componentes de la velocidad en el paso t+dt
            v_rot_new = p.omega * SVector(-r_c_wall[2], r_c_wall[1])
            v_cont_new = p.v + v_rot_new
            v_n_new_mag = dot(v_cont_new, n_wall)
            v_t_new = v_cont_new - v_n_new_mag * n_wall
            v_t_new_mag = norm(v_t_new)

            # Término normal viscoso disipado
            potencia_disipada_pared += 0.5 * gamma_n_wall * (v_n_old_mag^2 + v_n_new_mag^2)

            #= Término tangencial por fricción
            # Se utiliza la fuerza tangencial corregida en t + dt
            F_n_mag = max(k_n_wall * delta - gamma_n_wall * v_n_media_mag, 0.0)
            F_t_prueba = gamma_t_wall * v_t_media_mag
            F_t_coulomb = mu_wall * F_n_mag
            F_t_mag = min(F_t_prueba, F_t_coulomb)

            #Termina tangencial disipado
            potencia_disipada_pared += F_t_mag * v_t_media_mag
            =#
        end
    end

    energia_Mecanica[5] += potencia_disipada_pared * dt
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

# -- Funcion para exportar datos extraidos de la simulacion -- #
function guardar_datos(archivo::String, tiempo::Float64, energia_Mecanica::Vector{T}) where {T}
    energia_Cinetica_t = energia_Mecanica[1]            # Energia cinetica traslacional
    energia_Cinetica_r = energia_Mecanica[2]            # Energia cinetica rotacional
    energia_Cinetica = energia_Mecanica[3]              # Energia cinetica total
    energia_Pot_pared = energia_Mecanica[4]             # Energia potencial pared-particula
    energia_Disip_pared = energia_Mecanica[5]           # Energia disipada por colisiones pared-particula
    energia_Pot_part = energia_Mecanica[6]              # Energia potencial particula-particula
    energia_Disip_part = energia_Mecanica[7]            # Energia disipada por colisiones particula-particula
    trabajo_Inercial = energia_Mecanica[8]              # Trabajo inercial de la excitacion orbital
    energia_Cinetica_relativa = energia_Mecanica[9]     # Energia cinetica relativa al marco no inercial
    
    # La energia mecanica presente en los cuerpos instantaneamente
    energia_Mecanica_Total = energia_Cinetica_relativa + energia_Pot_pared + energia_Pot_part
    
    # La energia disipada total acumulada
    energia_Disipada_Total = energia_Disip_pared + energia_Disip_part
    
    # Balance = (Energía Mecánica + Energía Disipada) - Trabajo Inyectado
    # Este valor debería ser cercano a 0 en todo momento si el integrador de Verlet es estable.
    balance = (energia_Mecanica_Total + energia_Disipada_Total) - trabajo_Inercial

    open(archivo, "a") do io
        println(io,"$tiempo,","$energia_Cinetica_relativa,","$energia_Pot_pared,","$energia_Disip_pared,","$energia_Pot_part,","$energia_Disip_part,","$energia_Mecanica_Total,","$trabajo_Inercial,","$balance")
    end
end

# -- Funcion principal para simular el sistema de particulas -- #
function simular_sistema()
    # -- Parametros del sistema -- #
    dimension_Sistema = 2
    numero_Particulas = 90
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
        k_n_wall = 1.0e5, #
        gamma_n_wall = 1.0e3, #

        # -- Tangenciales -- #
        gamma_t_wall = 0.0e3, #
        mu_wall = 0.0, #

        # -- Parametros de colisiones particula-particula -- #
        
        # -- Normales -- #
        k_n = 1.0e5, #
        gamma_n = 1.0e3, #

        # -- Tangenciales -- #
        gamma_t = 1.0e3, #
        mu = 0.5 #
    )

    # -- Inicializacion del sistema de particulas -- #
    #= 
    sistema = generar_sistema(numero_Particulas, dimension_Sistema, R = radio_Particula, m = masa_Particula)
    generar_configuracion!(sistema, radio_Recipiente)
    =#
    # Sistema de particulas de prueba para verificar funciones o estabilidad de la simulacion
    
    sistema = [Particle(@SVector[5.0, 0.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5, 0.0, 0.0, 0.0, 0.5*masa_Particula*radio_Particula^2)#=, 
    Particle(@SVector[-1.0, 0.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5, 0.0, 0.0, 0.0, 0.5*masa_Particula*radio_Particula^2)=#]
    
    # -- Parametros de la ejecucion -- #
    dt = parametros_Generales.dt                    # En segundos (s)
    tiempo_total = 50.0                             # En segundos (s)
    pasos = round(Int, tiempo_total / dt)
    frecuencia_guardado = 100                       # Cada 0.01 segundos (s) se guardan los datos de la simulacion

    # -- Archivos de salida -- #
    # Nombre de los archivos de salida para los datos de la simulacion y para los datos extraidos de la simulacion
    # MALDITA SEA, ACUERDATE DE CAMBIAR EL NOMBRE DEL ARCHIVO DE SALIDA, NO TE VUELVAS A EQUIVOCAR
    archivo_salida_1 = "resultados/giro_swirling_172.xyz"
    archivo_salida_2 = "datos/datos_swirling_172.csv"
    # Apertura y limpieza de los archivos de salida
    open(archivo_salida_1, "w") do io end # Limpiar archivo si existe
    open(archivo_salida_2, "w") do io end # Limpiar archivo si existe

    # -- Arreglo para registrar la energia del sistema -- # 
    energia_Mecanica = [
        0.0,    # [1] Energia cinetica traslacional
        0.0,    # [2] Energia cinetica rotacional
        0.0,    # [3] Energia cinetica total

        # Variables de energia relacionadas con la interaccion particula-pared
        0.0,    # [4] Energia potencial total particula-pared
        0.0,    # [5] Energia disipada total por friccion particula-pared

        # Variables de energia relacionadas con la interaccion particula-particula
        0.0,    # [6] Energia potencial total particula-particula
        0.0,    # [7] Energia disipada total por friccion particula-particula

        # Variables de energia relacionadas con el marco no inercial
        0.0,    # [8] Trabajo inercial acumulado
        0.0     # [9] Energia cinetica relativa
        ]

    println("Iniciando movimiento swirled...")
    
    # -- Bucle principal de la simulacion -- #
    for paso in 1:pasos
        # -- Integracion de las ecuaciones de movimiento mediante el metodo de Velocity Verlet -- #
        velocity_verlet_step!(sistema, dt, fuerza_total!, dt * paso, parametros_Generales, energia_Mecanica, radio_Recipiente)
        
        # -- Guardado de los datos de la simulacion y de los datos extraidos de la simulacion -- #
        if paso % frecuencia_guardado == 0
            guardar_frame_xyz(archivo_salida_1, sistema, paso * dt, radio_Recipiente)
            guardar_datos(archivo_salida_2, paso * dt, energia_Mecanica)
        end
    end
    
    println("¡Simulación terminada!\nArchivos generados: $archivo_salida_1 y $archivo_salida_2")
end

simular_sistema()