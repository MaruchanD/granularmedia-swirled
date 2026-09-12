using Plots
using CSV
using DataFrames

# 1. Leer los datos desde el archivo .csv
df = CSV.read("datos/datos_swirling_172.csv", DataFrame; header=false)

tiempo            = df[:, 1]
energia_cinetica = df[:, 2]
energia_potencial = df[:, 3]
energia_disipada = df[:, 4]
energia_mecanica_total = df[:, 7]
trabajo_inercial = df[:, 8]
balance = df[:, 9]

# 2. Configuración de escalas, límites y leyendas
tipo_escala_x    = :identity
tipo_escala_y    = :identity

limites_x        = (0, 50)          # ej. (0, 50) o :auto
limites_y_p1     = :auto          # límites para el primer gráfico
limites_y_p2     = :auto          # límites para el segundo gráfico

posicion_leyenda = :outertop

# 3. Primer gráfico (Energía Cinética y Potencial)
p1 = plot(tiempo, energia_cinetica,
          label="Energia Cinetica",
          ylabel="Energía (E)",
          title="Componentes de Energía",
          linewidth=2,
          xscale=tipo_escala_x,
          yscale=tipo_escala_y,
          xlims=limites_x,
          ylims=limites_y_p1,
          legend=posicion_leyenda
)
#=
plot!(p1, tiempo, energia_potencial,
      label="Energía Potencial",
      linewidth=2
)

plot!(p1, tiempo, energia_mecanica_total,
      label="Energía Mecanica Total",
      linewidth=2
)
=#
plot!(p1, tiempo, energia_disipada,
      label="Energia Disipada",
      linewidth=2
)

plot!(p1, tiempo, trabajo_inercial,
      label="Trabajo Inercial",
      linewidth=2
)

# 4. Segundo gráfico (Energía Total)
p2 = plot(tiempo, balance,
          label="Balance Energetico Total",
          xlabel="Tiempo (t)",
          ylabel="Energía (E)",
          title="Comparacion",
          linewidth=2,
          xscale=tipo_escala_x,
          yscale=tipo_escala_y,
          xlims=limites_x,
          ylims=limites_y_p2,
          legend=posicion_leyenda)

# 5. Combinar ambos paneles en una sola imagen (2 filas, 1 columna)
# link=:x sincroniza el rango del eje temporal en ambas gráficas
grafico_final = plot(p1, p2, 
                     layout=(2, 1), 
                     link=:x, 
                     size=(800, 600))

# 6. Guardar la imagen combinada
savefig(grafico_final, "graficos/balance_y_energia_cinetica_1-particula.png")