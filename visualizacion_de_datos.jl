using Plots
using CSV
using DataFrames

# 1. Leer los datos desde el archivo .csv
# Si el archivo NO tiene encabezados (nombres de columnas en la primera fila), usa header=false:
df = CSV.read("datos/datos_swirling_170.csv", DataFrame; header=false)

# 2. Extraer las columnas por posición (o por nombre si tuviera encabezados, ej: df.tiempo)
tiempo            = df[:, 1]
energia_cinetica_traslacional = df[:, 2]
energia_cinetica_rotacional = df[:, 3]
energia_potencial_particula_pared = df[:, 4]
energia_disipada_particula_pared = df[:, 5]
energia_potencial_particula_particula = df[:, 6]
energia_disipada_particula_particula = df[:, 7]
energia_total = df[:, 8]

# Posición de la leyenda:
# Ejemplos dentro: :topright, :topleft, :bottomleft, :bottomright, :best
# Ejemplos fuera: :outertopright, :outerright, :outertop
posicion_leyenda = :outertop

# 3. Configuración de escalas y límites
tipo_escala_x = :identity   # :identity (lineal), :log10, :ln
tipo_escala_y = :identity

limites_x = (4, 10)         # Define (xmin, xmax) o usa :auto
limites_y = :auto #(-1, 10)      # Define (ymin, ymax) o usa :auto

# 4. Graficar
plot(tiempo, energia_cinetica_traslacional,
     label="Energía Cinetica Traslacional",
     xlabel="Tiempo (t)",
     ylabel="Energía (E)",
     title="Evolución de las Energías en el Tiempo",
     linewidth=2,
     xscale=tipo_escala_x,
     yscale=tipo_escala_y,
     xlims=limites_x,
     ylims=limites_y,
     legend=posicion_leyenda  # <--- Control de la posición
)

#plot!(tiempo, energia_cinetica_rotacional, label="Energía Cinetica Rotacional", linewidth=2)

#plot!(tiempo, energia_potencial_particula_pared, label="Energía Potencial (Partícula-Pared)", linewidth=2)

#plot!(tiempo, energia_potencial_particula_particula, label="Energía Potencial (Partícula-Partícula)", linewidth=2)

#plot!(tiempo, energia_disipada_particula_pared, label="Energía Disipada (Partícula-Pared)", linewidth=2)

#plot!(tiempo, energia_disipada_particula_particula, label="Energía Disipada (Partícula-Partícula)", linewidth=2)

#plot!(tiempo, energia_total, label="Energía Total", linewidth=2)

# 5. Guardar el gráfico
savefig("graficos/evolucion_energia_cinetica_traslacional_1_particula_sin_roce_tangencial.png")