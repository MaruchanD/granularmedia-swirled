# ==============================================================================
# Script de automatización para barrido de simulaciones por número de partículas
# ==============================================================================

# 1. Incluir el archivo base que contiene las funciones y el modelo
# (Asegúrate de que el archivo con tu código se llame así o ajusta la ruta)
include("simulacion_swirling_v0.200.jl")

function ejecutar_barrido()
    println("=====================================================")
    println("   BARRIDO DE SIMULACIONES: MATERIA GRANULAR (DEM)   ")
    println("=====================================================")

    # -- Entrada interactiva de los límites del rango -- #
    print("Introduce el número de partículas INICIAL (N_inicio): ")
    n_inicio = parse(Int, readline())

    print("Introduce el número de partículas FINAL (N_fin): ")
    n_fin = parse(Int, readline())

    # Validación básica del rango
    if n_inicio > n_fin
        println("\n[Error] El número inicial ($n_inicio) no puede ser mayor que el final ($n_fin).")
        return
    end

    if n_inicio <= 0
        println("\n[Error] El número de partículas debe ser un entero positivo.")
        return
    end

    # Crear carpetas de destino automáticamente si no existen
    isdir("resultados") || mkpath("resultados")
    isdir("datos") || mkpath("datos")

    total_sims = (n_fin - n_inicio) + 1
    println("\n-> Se ejecutarán $total_sims simulación(es) en total (desde N = $n_inicio hasta N = $n_fin).")
    println("-----------------------------------------------------")

    tiempo_inicio_total = time()

    # -- Bucle de ejecución para cada valor de N (de uno en uno) -- #
    for (idx, n_particulas) in enumerate(n_inicio:n_fin)
        println("\n[$idx/$total_sims] Ejecutando simulación para N = $n_particulas...")
        t_iter_inicio = time()

        # Llamar a la función simular_sistema pasándole explícitamente el valor de N
        simular_sistema(n_particulas)

        t_iter_fin = time()
        duracion_iter = round(t_iter_fin - t_iter_inicio, digits=2)
        println("-> Simulación para N = $n_particulas completada en $(duracion_iter) s.")
    end

    tiempo_total = round((time() - tiempo_inicio_total) / 60.0, digits=2)
    println("\n=====================================================")
    println(" ¡Barrido completo finalizado con éxito!")
    println(" Tiempo total acumulado: $tiempo_total minutos.")
    println("=====================================================")
end

# Ejecutar el barrido al correr el script
ejecutar_barrido()