# ==============================================================================
# Script para comparar curvas f vs N para distintas amplitudes (A)
# ==============================================================================

using DelimitedFiles
using Plots

function graficar_comparacion_amplitudes(directorio_datos::String = "datos")
    if !isdir(directorio_datos)
        error("El directorio '$directorio_datos' no existe.")
    end

    # Patrón para identificar los archivos de resumen y capturar amplitud (A) y frecuencia (F)
    patron_resumen = r"^resumen_macroscopico_A([\d\.]+)_F([\d\.]+)\.csv$"
    archivos = readdir(directorio_datos)

    # Estructuras para almacenar los datos extraídos
    datos_por_amplitud = Dict{Float64, Tuple{Vector{Int}, Vector{Float64}, String}}()

    println("Leyendo archivos de resumen en '$directorio_datos'...")

    for nombre_archivo in archivos
        coincidencia = match(patron_resumen, nombre_archivo)
        
        if coincidencia !== nothing
            amp_val = parse(Float64, coincidencia.captures[1])
            frec_val_str = coincidencia.captures[2]
            ruta = joinpath(directorio_datos, nombre_archivo)

            # header=true omite la primera fila ("N,frecuencia_promedio_Hz,inercia_promedio")
            matriz, encabezados = readdlm(ruta, ',', Float64; header=true)

            if size(matriz, 1) > 0 && size(matriz, 2) >= 2
                N_vals = round.(Int, matriz[:, 1])
                f_vals = matriz[:, 2]

                # Ordenar por N de forma ascendente
                orden = sortperm(N_vals)
                datos_por_amplitud[amp_val] = (N_vals[orden], f_vals[orden], frec_val_str)
                println("-> Cargado: $nombre_archivo (A = $amp_val, F = $frec_val_str Hz)")
            else
                @warn "El archivo $nombre_archivo no contiene datos válidos."
            end
        end
    end

    if isempty(datos_por_amplitud)
        println("No se encontraron archivos coincidentes con el patrón esperado.")
        return
    end

    # Ordenar las curvas por valor creciente de amplitud para una leyenda limpia
    amplitudes_ordenadas = sort(collect(keys(datos_por_amplitud)))

    # Inicializar el lienzo del gráfico
    p = plot(
        title = "Frecuencia de rotación vs Número de partículas",
        xlabel = "Número de partículas (N)",
        ylabel = "Frecuencia ⟨f⟩ [Hz]",
        grid = true,
        legend = :topright,
        size = (850, 550)
    )

    # Línea de referencia en f = 0 (cruce rotación a reptación)
    hline!(p, [0.0], linestyle = :dash, color = :black, label = "f = 0 (Transición)")

    # Trazar cada serie correspondiente a una amplitud distinta
    for amp in amplitudes_ordenadas
        N_vals, f_vals, _ = datos_por_amplitud[amp]
        plot!(
            p,
            N_vals,
            f_vals,
            seriestype = :scatter,
            line = true,
            marker = :circle,
            markersize = 4,
            label = "A = $amp cm"
        )
    end

    display(p)

    nombre_salida = "comparacion_frecuencia_vs_N_amplitudes.png"
    savefig(p, nombre_salida)
    println("\n-> Gráfico comparativo guardado exitosamente como '$nombre_salida'.")
end

graficar_comparacion_amplitudes()