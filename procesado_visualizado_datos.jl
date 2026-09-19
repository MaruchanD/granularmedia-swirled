# ==============================================================================
# Script para análisis y graficado macroscópico en función de N
# ==============================================================================

using DelimitedFiles
using Plots

function procesar_y_graficar_datos(directorio_datos::String = "datos")
    if !isdir(directorio_datos)
        error("El directorio '$directorio_datos' no existe.")
    end

    # Expresión regular para filtrar y validar el formato de los archivos:
    # Coincide con: macroscopicos_N<int>_A<num>_F<num>.csv o datos_macroscopicos_N...
    patron_archivo = r"^(?:datos_)?macroscopicos_N(\d+)_A([\d\.]+)_F([\d\.]+)\.csv$"

    lista_archivos = readdir(directorio_datos)
    
    # Arreglos para almacenar los resultados consolidados
    valores_N = Int[]
    frecuencia_promedio = Float64[]
    inercia_promedio = Float64[]

    # Variables para almacenar la amplitud y frecuencia leídas del nombre
    amplitud_str = ""
    frecuencia_str = ""

    println("Buscando y procesando archivos en '$directorio_datos'...")

    for nombre_archivo in lista_archivos
        coincidencia = match(patron_archivo, nombre_archivo)
        
        if coincidencia !== nothing
            # Extraer amplitud y frecuencia del primer archivo que coincida
            if isempty(amplitud_str)
                amplitud_str = coincidencia.captures[2]
                frecuencia_str = coincidencia.captures[3]
            end

            ruta_completa = joinpath(directorio_datos, nombre_archivo)
            
            # Leer el archivo delimitado por comas
            datos = readdlm(ruta_completa, ',', Float64)

            # Validar que el archivo contenga datos y al menos 4 columnas
            if size(datos, 1) > 0 && size(datos, 2) >= 4
                # Columna 2: Número de partículas (tomamos la primera fila como referencia)
                N_particulas = round(Int, datos[1, 2])
                
                # Columna 3: Frecuencia de rotación (Hz)
                # Columna 4: Momento de inercia
                # Promediamos a lo largo de las filas (dimensión temporal)
                f_mean = sum(datos[:, 3]) / size(datos, 1)
                I_mean = sum(datos[:, 4]) / size(datos, 1)

                push!(valores_N, N_particulas)
                push!(frecuencia_promedio, f_mean)
                push!(inercia_promedio, I_mean)
                
                println("-> Procesado: $nombre_archivo (N = $N_particulas)")
            else
                @warn "El archivo $nombre_archivo está vacío o no contiene las columnas requeridas."
            end
        end
    end

    if isempty(valores_N)
        println("No se encontraron archivos válidos con el patrón indicado.")
        return
    end

    # Ordenar los datos en función de N de menor a mayor
    orden = sortperm(valores_N)
    valores_N = valores_N[orden]
    frecuencia_promedio = frecuencia_promedio[orden]
    inercia_promedio = inercia_promedio[orden]

    # =========================================================================================================
    # Exportacion de los datos obtenidos como funcion del numero de particulas, la frecuencia  y la amplitud
    # =========================================================================================================
    archivo_resumen_csv = joinpath(directorio_datos, "resumen_macroscopico_A$(amplitud_str)_F$(frecuencia_str).csv")
    
    open(archivo_resumen_csv, "w") do io
        # Encabezado descriptivo
        println(io, "N,frecuencia_promedio_Hz,inercia_promedio")
        # Escritura de los datos ordenados
        for i in 1:length(valores_N)
            println(io, "$(valores_N[i]),$(frecuencia_promedio[i]),$(inercia_promedio[i])")
        end
    end
    println("\n-> Datos consolidados guardados en: '$archivo_resumen_csv'")

    # ==============================================================================
    # Generación de gráficos (subplots compartiendo el eje X)
    # ==============================================================================
    
    # Gráfico superior: Frecuencia de rotación promedio vs N
    p1 = plot(
        valores_N, 
        frecuencia_promedio,
        seriestype = :scatter,
        line = true,
        marker = :circle,
        color = :blue,
        label = "⟨f⟩ (Hz)",
        ylabel = "Frecuencia [Hz]",
        title = "Dinámica macroscópica del clúster",
        grid = true,
        legend = :topright
    )
    # Línea horizontal en f = 0 para evidenciar la transición rotación / reptación
    hline!(p1, [0.0], linestyle = :dash, color = :black, label = "f = 0")

    # Gráfico inferior: Momento de inercia promedio vs N
    p2 = plot(
        valores_N, 
        inercia_promedio,
        seriestype = :scatter,
        line = true,
        marker = :square,
        color = :red,
        label = "⟨I⟩",
        xlabel = "Número de partículas (N)",
        ylabel = "Momento de Inercia",
        grid = true,
        legend = :topleft
    )

    # Combinación de ambos gráficos compartiendo eje horizontal y alineados verticalmente
    figura_final = plot(
        p1, 
        p2, 
        layout = (2, 1), 
        link = :x,            # Enlaza el rango del eje x entre ambos subplots
        size = (800, 700)
    )

    # Mostrar en pantalla y guardar en disco
    display(figura_final)
    savefig(figura_final, "resultados_macroscopicos_vs_N_(A$(amplitud_str)_F$(frecuencia_str)).png")
    println("\nGráfica guardada exitosamente como 'resultados_macroscopicos_vs_N_(A$(amplitud_str)_F$(frecuencia_str)).png'.")
end

# Ejecución
procesar_y_graficar_datos()