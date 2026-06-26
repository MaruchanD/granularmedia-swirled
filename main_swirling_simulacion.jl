using StaticArrays
using LinearAlgebra

# ==========================================
#   1. DEFINICIÓN DE LA ESTRUCTURA
# ==========================================

mutable struct Particula{N, T}
    r::SVector{N, T}        #Posicion
    v::SVector{N, T}        #Velocidad
    a::SVector{N, T}        #Aceleracion
    m::T                    #Masa
    radio::T                #Radio
end

# =================================================
#   2. ALGORITMO DE INTEGRACIÓN (VELOCITY-VERLET)
# =================================================

function velocity_verlet_step!(particulas::Vector{Particula{N, T}}, dt::T, calc_fuerzas!) where {N, T}

    # PASO 1 y 2: Se actualizan la posicion y la velocidad a medio paso.

    for p in particulas
        p.r = p.r + p.v * dt + 0.5 * p.a * dt^2
        p.v = p.v + 0.5 * p.a * dt
    end

    # PASO 3: Se calculan las nuevas fuerzas/aceleraciones.
    # Se pasan los arreglos mutados a la funcion de fuerzas

    calc_fuerzas!(particulas)

    # PASO 4: Completamos la actualizacion del paso de velocidad

    for p in particulas
        p.v = p.v + 0.5 * p.a * dt
    end
end

# =================================================
# 3. FUNCIÓN DE FUERZAS (Linear spring-dashpot model)
# =================================================

function spring_dashpot!(particulas::Vector{Particula{N, T}}) where {N, T}
    # Constantes del modelo o Parametros del material
    k = 1000.0  # Constante de resorte (N/m), dependiendo de la dureza del material
    gamma = 10.0  # Coeficiente de amortiguamiento (Ns/m) dependiendo de la viscosidad del material

    # PASO 1. Iniciar las aceleraciones a cero
    for p in particulas
        p.a = @SVector zeros(T, N)
    end

    num_p = length(particulas)

    # PASO 2. Bucle sobre pares unicos de particulas O(N^2)
    for i in 1:num_p
        p_i = particulas[i]

        for j in (i+1):num_p
            p_j = particulas[j]

            # Vector relativo de posiciones (de j a i)
            r_ij = p_j.r - p_i.r

            # Magnitud del vector r_ij (o la distancia d entre los centros de las particulas)
            d = norm(r_ij)

            # Suma de los radios de las particulas
            suma_radios = p_i.radio + p_j.radio

            # Verificacion de solapamiento (colision)
            if d < suma_radios && d > 0.0
                delta = suma_radios - d     # Se define el solapamiento delta
                n_ij = r_ij / d             # Vector normal unitario

                # Velocidad relativa
                v_ij = p_j.v - p_i.v
                vn = dot(v_ij, n_ij)        # Proyeccion normal de la velocidad relativa

                #Calculo de la magnitud de la fuerza normal (Resorte-Amortiguador)
                Fn_mag = k * delta - gamma * vn
                Fn_mag = max(Fn_mag, 0.0)   #Evitar fuerzas atractivas artificiales, debido a que la disipacion se puede hacer grande

                # Vector de fuerza normal
                Fn_vec = Fn_mag * n_ij

                # Actualizar aceleraciones de las particulas (F = m * a -> a = F / m)
                # Se acumula fuerza directamente sobre el SVector mediante la 3era ley de Newton
                particulas[i].a += Fn_vec / p_i.m
                particulas[j].a -= Fn_vec / p_j.m
            end
        end
    end
end

function aplicar_fuerzas_contenedor_circular!(particulas::Vector{Particula{N, T}}, r_contenedor::T) where {N, T}
    k_pared = 1.0e5
    gamma_pared = 1.0e2

    for p in particulas
        # 1. Distancia desde el centro (0.0, 0.0)
        d = norm(p.r)
        # 2. Calcular el solapamiento con la pared circular
        delta = d + p.radio - r_contenedor
        # 3. Evaluar si hay solapamiento
        if delta > 0.0
            # 4. Calcular la fuerza normal de la pared
            n_pared = -p.r / d              # Vector normal unitario hacia el interior del recipiente
            v_n = dot(p.v, n_pared)         # Componente normal de la velocidad

            # Calcular la magnitud de la fuerza normal (Resorte-Amortiguador)
            # Se evalua con max() para evitar la aparicion de fuerzas de atraccion.
            Fn_mag = max(k_pared * delta - gamma_pared * v_n, 0.0)
            Fn_vec = Fn_mag * n_pared

            # 5. Actualizar aceleración de la partícula
            p.a += Fn_vec / p.m
        end
    end
end

