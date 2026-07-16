# Intento para agregar fuerzas tangenciales entre particulas. #
# Se usa un enfoque viscoso regularizado. Se modela el contacto como un amortiguador viscoso, pero truncado
# por el limite de Coulomb. De aqui el calculo Ft_mag= min(Ft_prueba, F_limite_coulomb). Esto con el objeto
# de no evaluar esta fuerza de forma instantanea y tambien no modificar todo la estructura de datos utilizada
# hasta ahora.
using StaticArrays
using LinearAlgebra

mutable struct Particle{N, T}
    r::SVector{N, T}      # Posición
    v::SVector{N, T}      # Velocidad
    a::SVector{N, T}      # Aceleración
    mass::T               # Masa
    radius::T             # Radio
end

function velocity_verlet_step!(particles::Vector{Particle{N, T}}, dt::T, calc_forces!, tiempo::T) where {N, T}
    for p in particles
        p.r = p.r + p.v * dt + 0.5 * p.a * dt^2
        p.v = p.v + 0.5 * p.a * dt
    end
    
    calc_forces!(particles, tiempo + dt)
    
    for p in particles
        p.v = p.v + 0.5 * p.a * dt
    end
end

function excitacion_orbital_rampa!(particles::Vector{Particle{N,T}}, tiempo::T) where {N, T}
    # Parametros para la excitacion
    amplitud = 3.0                      # Radio de la excitacion (1.0 a 5.0 cm)
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

function contenedor_circular!(particles::Vector{Particle{N, T}}) where {N, T}
    # Parametros de la interaccion con el contenedor
    radio_contenedor = 5.0
    # Parametros normales
    kn_wall = 1.0e5
    gamma_wall = 5.0e3
    # Parametros tangenciales
    mu_pared = 0.3 
    gamma_t_pared = 8.0e2

    for p in particles
        d = norm(p.r)
        delta = d + p.radius - radio_contenedor

        # Aplicacion de las fuerzas de contacto
        if delta > 0.0
            n_wall = -p.r / d               # Direccion normal al centro
            v = p.v                         # Velocidad de la particula
            
            v_n_mag = dot(p.v, n_wall)      # Magnitud de la velocidad en la direccion normal
            v_n_vec = v_n_mag * n_wall      # Componente normal de la velocidad
            
            v_t_vec = v - v_n_vec           # Componente tangencial de la velocidad
            v_t_mag = norm(v_t_vec)         # Magnitud de la velocidad en la direccion tangencial

            # Fuerza de contacto normal con la pared
            F_n_mag = max(kn_wall * delta - gamma_wall * v_n_mag, 0.0)
            F_n_vec = F_n_mag * n_wall
            # Aplicacion de la fuerza normal
            p.a += F_n_vec / p.mass
            
            #Fuerza de contacto tangencial con la pared
            if v_t_mag > 0.0
                t_wall = v_t_vec / v_t_mag

                F_tangencial_prueba = gamma_t_pared * v_t_mag
                F_limite_coulomb = mu_pared * F_n_mag

                F_t_mag = min(F_tangencial_prueba, F_limite_coulomb)
                F_t_vec = -F_t_mag * t_wall
                # Aplicacion de la fuerza tangencial
                p.a += F_t_vec / p.mass
            end
        end
    end
end


function contacto_particulas!(particles::Vector{Particle{N, T}}) where {N, T}
    # Parametros normales
    k_n = 1.0e5
    gamma_n = 5.0e1

    # Parametros tangenciales
    mu = 0.5
    gamma_t = 1.0e3

    num_p = length(particles)

    for i in 1:num_p
        p_i = particles[i]
        for j in (1+i):num_p
            p_j = particles[j]

            r_ij = p_i.r - p_j.r
            d = norm(r_ij)
            suma_radios = p_i.radius + p_j.radius

            if d < suma_radios && d > 0.0
                # Vector Normal Unitario
                n_ij = r_ij / d

                # Descomposicion de la velocidad relativa en componentes normal y tangencial
                v_ij = p_i.v - p_j.v                # Velocidad relativa
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

                #Fuerza Tangencial
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

                    #Aplicacion de la fuerza tangencial
                    particles[i].a += Ft_vec / p_i.mass
                    particles[j].a -= Ft_vec / p_j.mass
                end
            end
        end
    end
end

function fuerza_total!(particles::Vector{Particle{N, T}}, tiempo::T) where {N, T}
    # Reiniciar aceleraciones
    for p in particles
        p.a = @SVector zeros(T, N)
    end

    # Aplicar todas las fuerzas
    excitacion_orbital_rampa!(particles, tiempo)
    contenedor_circular!(particles)
    contacto_particulas!(particles)
end

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

function simular_sistema()
    sistema = [
        Particle(@SVector[0.0, 0.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5),
        Particle(@SVector[2.0, 1.0], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5),
        Particle(@SVector[-2.5, -1.5], @SVector[0.0, 0.0], @SVector[0.0, 0.0], 1.0, 0.5)
    ]
    # Parametros
    radio_contenedor = 5.0

    dt = 1.0e-4
    tiempo_total = 10.0 
    pasos = round(Int, tiempo_total / dt)
    frecuencia_guardado = 100 
    
    archivo_salida = "resultados/giro_swirling_130.xyz"
    open(archivo_salida, "w") do io end # Limpiar archivo si existe

    println("Iniciando movimiento swirled...")
    
    for paso in 1:pasos
        velocity_verlet_step!(sistema, dt, fuerza_total!, dt * paso)
        
        if paso % frecuencia_guardado == 0
            guardar_frame_xyz!(archivo_salida, sistema, paso * dt, radio_contenedor)
        end
    end
    
    println("¡Simulación terminada! Archivo generado: $archivo_salida")
end

simular_sistema()