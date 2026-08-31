using LinearAlgebra
using Plots

"""
Genera una configuración densa sin superposiciones utilizando 
un algoritmo de relajación (dinámica sobreamortiguada).
"""
function generar_configuracion_densa(N::Int, R_container::Float64, r_particle::Float64; max_iters::Int=20000, dt::Float64=0.1, tol::Float64=1e-6)
    
    # Validación física del empaquetamiento (Límite en 2D es aprox 0.84 - Random Close Packing)
    area_particulas = N * pi * r_particle^2
    area_recipiente = pi * R_container^2
    fraccion_empaquetamiento = area_particulas / area_recipiente
    
    println("Fracción de empaquetamiento requerida: ", round(fraccion_empaquetamiento, digits=3))
    if fraccion_empaquetamiento > 0.82
        @warn "La fracción de empaquetamiento está muy cerca o supera el límite teórico de Random Close Packing en 2D. Es posible que el sistema no logre converger."
    end

    # 1. Inserción inicial aleatoria (se permiten superposiciones)
    pos = zeros(Float64, N, 2)
    R_eff = R_container - r_particle
    
    for i in 1:N
        r = R_eff * sqrt(rand())
        theta = 2 * pi * rand()
        pos[i, 1] = r * cos(theta)
        pos[i, 2] = r * sin(theta)
    end

    # 2. Bucle de Relajación
    diameter = 2.0 * r_particle
    diam_sq = diameter^2
    displacements = zeros(Float64, N, 2)
    
    for iter in 1:max_iters
        max_overlap = 0.0
        fill!(displacements, 0.0) # Reiniciar desplazamientos
        
        # Interacción Partícula-Partícula
        for i in 1:(N-1)
            for j in (i+1):N
                dx = pos[i, 1] - pos[j, 1]
                dy = pos[i, 2] - pos[j, 2]
                dist_sq = dx^2 + dy^2
                
                if dist_sq < diam_sq && dist_sq > 0.0
                    dist = sqrt(dist_sq)
                    overlap = diameter - dist
                    max_overlap = max(max_overlap, overlap)
                    
                    # Fuerza repulsiva en dirección normal
                    fx = (dx / dist) * overlap * 0.5
                    fy = (dy / dist) * overlap * 0.5
                    
                    displacements[i, 1] += fx
                    displacements[i, 2] += fy
                    displacements[j, 1] -= fx
                    displacements[j, 2] -= fy
                end
            end
        end
        
        # Interacción Partícula-Pared
        for i in 1:N
            r_sq = pos[i, 1]^2 + pos[i, 2]^2
            if r_sq > R_eff^2
                dist = sqrt(r_sq)
                overlap = dist - R_eff
                max_overlap = max(max_overlap, overlap)
                
                # Fuerza de restauración hacia el centro
                fx = -(pos[i, 1] / dist) * overlap
                fy = -(pos[i, 2] / dist) * overlap
                
                displacements[i, 1] += fx
                displacements[i, 2] += fy
            end
        end
        
        # Actualizar posiciones
        pos .+= dt .* displacements
        
        # Criterio de convergencia: si la superposición máxima es imperceptible
        if max_overlap < tol
            println("Convergencia alcanzada en la iteración $iter.")
            return pos
        end
    end
    
    @warn "El sistema no convergió completamente tras $max_iters iteraciones. Superposición máxima residual: $max_overlap"
    return pos
end

# ==========================================
# Ejecución y Visualización
# ==========================================

N_particulas = 40
Radio_recipiente = 5.0
Radio_particula = 0.5

# Generar posiciones
coords = generar_configuracion_densa(N_particulas, Radio_recipiente, Radio_particula)

# Graficar el resultado
gr() # Usar el backend estándar
plt = plot(aspect_ratio=:equal, legend=false, title="Estado Inicial Maduro", grid=true)

# Dibujar el recipiente (aproximado con muchos puntos)
theta = range(0, 2pi, length=200)
plot!(plt, Radio_recipiente .* cos.(theta), Radio_recipiente .* sin.(theta), color=:black, linewidth=2)

# Dibujar las partículas
for i in 1:N_particulas
    x, y = coords[i, 1], coords[i, 2]
    # plot! en Julia requiere pasar los puntos del contorno para dibujar círculos eficientemente
    plot!(plt, x .+ Radio_particula .* cos.(theta), y .+ Radio_particula .* sin.(theta), fill=(0, :blue), fillalpha=0.6, linecolor=:black)
end

display(plt)