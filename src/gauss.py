import os

carpeta_proyecto = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def cargar_datos():
    ruta_matriz = os.path.join(carpeta_proyecto, "datos", "matriz_a.csv")
    ruta_vector = os.path.join(carpeta_proyecto, "datos", "vector_b.csv")

    matriz = []
    archivo = open(ruta_matriz)
    for linea in archivo:
        numeritos = linea.strip().split(",")
        filita = []
        for numerito in numeritos:
            filita.append(float(numerito))
        matriz.append(filita)
    archivo.close()

    vector = []
    archivo = open(ruta_vector)
    for linea in archivo:
        vector.append(float(linea.strip()))
    archivo.close()

    return matriz, vector

def mayor_valor_absoluto(matriz):
    tamano = len(matriz)
    mayor = 0.0
    for fila in range(tamano):
        for columna in range(tamano):
            if abs(matriz[fila][columna]) > mayor:
                mayor = abs(matriz[fila][columna])
    return mayor

def calcular_escalas(matriz):
    # escala de cada fila = el numero mas grande (en valor absoluto) de esa fila
    tamano = len(matriz)
    escalas = []
    comparaciones = 0
    for fila in range(tamano):
        mayor_de_la_fila = abs(matriz[fila][0])
        for columna in range(1, tamano):
            comparaciones = comparaciones + 1
            if abs(matriz[fila][columna]) > mayor_de_la_fila:
                mayor_de_la_fila = abs(matriz[fila][columna])
        escalas.append(mayor_de_la_fila)
    return escalas, comparaciones

def elegir_fila_pivote(matriz_copia, escalas, fila_pivote):
    # busca la fila (desde fila_pivote hacia abajo) con mayor |a(i,k)| / escala(i)
    tamano = len(matriz_copia)
    mejor_fila = fila_pivote
    mejor_valor = abs(matriz_copia[fila_pivote][fila_pivote]) / escalas[fila_pivote]
    comparaciones = 0
    divisiones = 1
    for fila_candidata in range(fila_pivote + 1, tamano):
        valor = abs(matriz_copia[fila_candidata][fila_pivote]) / escalas[fila_candidata]
        divisiones = divisiones + 1
        comparaciones = comparaciones + 1
        if valor > mejor_valor:
            mejor_valor = valor
            mejor_fila = fila_candidata
    return mejor_fila, comparaciones, divisiones

def eliminacion_gaussiana(matriz, vector, con_pivoteo):
    tamano = len(matriz)

    matriz_copia = []       # copio la matriz parte por parte
    for fila in range(tamano):
        filita_copia = []
        for columna in range(tamano):
            filita_copia.append(matriz[fila][columna])
        matriz_copia.append(filita_copia)

    vector_copia = []       # copio el vector numero por numero
    for fila in range(tamano):
        vector_copia.append(vector[fila])
    multiplicadores = []    # aqui guardo cada multiplicador (la matriz L)
    for fila in range(tamano):
        filita = []
        for columna in range(tamano):
            filita.append(0.0)
        multiplicadores.append(filita)
        
    orden_filas = []        # que fila original quedo en cada posicion
    for fila in range(tamano):
        orden_filas.append(fila)

    multdiv_por_paso = []       # contador de multiplicaciones y divisiones
    sumres_por_paso = []        # contador de sumas y restas
    comparaciones_por_paso = []     # busqueda del pivote (sobrecosto)
    divisiones_escala_por_paso = []     # |a(i,k)|/escala (sobrecosto)
    intercambios_por_paso = []      # cambios de fila (sobrecosto)
    for paso in range(tamano - 1):
        multdiv_por_paso.append(0)
        sumres_por_paso.append(0)
        comparaciones_por_paso.append(0)
        divisiones_escala_por_paso.append(0)
        intercambios_por_paso.append(0)

    comparaciones_escalas = 0
    if con_pivoteo:
        escalas, comparaciones_escalas = calcular_escalas(matriz_copia)

    mayor_original = mayor_valor_absoluto(matriz_copia)
    mayor_visto = mayor_original        # para el factor de crecimiento
    multiplicador_maximo = 0.0
    traza = []              # aqui guardo la matriz despues de cada paso

    for fila_pivote in range(tamano - 1):       # cada vueltita usa una fila como pivote

        if con_pivoteo:
            mejor_fila, comparaciones, divisiones = elegir_fila_pivote(matriz_copia, escalas, fila_pivote)
            comparaciones_por_paso[fila_pivote] = comparaciones
            divisiones_escala_por_paso[fila_pivote] = divisiones

            if mejor_fila != fila_pivote:       # intercambio las filas
                matriz_copia[fila_pivote], matriz_copia[mejor_fila] = matriz_copia[mejor_fila], matriz_copia[fila_pivote]
                vector_copia[fila_pivote], vector_copia[mejor_fila] = vector_copia[mejor_fila], vector_copia[fila_pivote]
                multiplicadores[fila_pivote], multiplicadores[mejor_fila] = multiplicadores[mejor_fila], multiplicadores[fila_pivote]
                escalas[fila_pivote], escalas[mejor_fila] = escalas[mejor_fila], escalas[fila_pivote]
                orden_filas[fila_pivote], orden_filas[mejor_fila] = orden_filas[mejor_fila], orden_filas[fila_pivote]
                intercambios_por_paso[fila_pivote] = 1

        pivote = matriz_copia[fila_pivote][fila_pivote]

        for fila_de_abajo in range(fila_pivote + 1, tamano):

            numero_a_eliminar = matriz_copia[fila_de_abajo][fila_pivote]
            multiplicador = numero_a_eliminar / pivote
            multdiv_por_paso[fila_pivote] = multdiv_por_paso[fila_pivote] + 1      # conte 1 division
            multiplicadores[fila_de_abajo][fila_pivote] = multiplicador     # lo guardo para reutilizarlo

            if abs(multiplicador) > multiplicador_maximo:
                multiplicador_maximo = abs(multiplicador)

            matriz_copia[fila_de_abajo][fila_pivote] = 0.0      #este numero ya se queda en cero, lo pongo directo

            for columna in range(fila_pivote + 1, tamano):      #esto actualiza el resto de la fila
                numero_de_arriba = matriz_copia[fila_pivote][columna]
                matriz_copia[fila_de_abajo][columna] = matriz_copia[fila_de_abajo][columna] - multiplicador * numero_de_arriba
                multdiv_por_paso[fila_pivote] = multdiv_por_paso[fila_pivote] + 1   # 1 multiplicacion
                sumres_por_paso[fila_pivote] = sumres_por_paso[fila_pivote] + 1     # 1 resta

                if abs(matriz_copia[fila_de_abajo][columna]) > mayor_visto:
                    mayor_visto = abs(matriz_copia[fila_de_abajo][columna])

            vector_copia[fila_de_abajo] = vector_copia[fila_de_abajo] - multiplicador * vector_copia[fila_pivote]       #al vector le pasa exactamente lo mismo, con el mismo multiplicador
            multdiv_por_paso[fila_pivote] = multdiv_por_paso[fila_pivote] + 1       # 1 multiplicacion
            sumres_por_paso[fila_pivote] = sumres_por_paso[fila_pivote] + 1         # 1 resta


        # foto de la matriz y el vector despues de este paso
        foto_matriz = []
        for fila in range(tamano):
            filita = []
            for columna in range(tamano):
                filita.append(matriz_copia[fila][columna])
            foto_matriz.append(filita)
        foto_vector = []
        for fila in range(tamano):
            foto_vector.append(vector_copia[fila])
        traza.append([foto_matriz, foto_vector, orden_filas[fila_pivote] + 1])
    for fila in range(tamano):      # la diagonal de L son unos
        multiplicadores[fila][fila] = 1.0
    resultado = {}
    resultado["matriz_triangular"] = matriz_copia
    resultado["vector_nuevo"] = vector_copia
    resultado["orden_filas"] = orden_filas
    resultado["multdiv_por_paso"] = multdiv_por_paso
    resultado["sumres_por_paso"] = sumres_por_paso
    resultado["comparaciones_por_paso"] = comparaciones_por_paso
    resultado["divisiones_escala_por_paso"] = divisiones_escala_por_paso
    resultado["intercambios_por_paso"] = intercambios_por_paso
    resultado["comparaciones_escalas"] = comparaciones_escalas
    resultado["multiplicador_maximo"] = multiplicador_maximo
    resultado["factor_crecimiento"] = mayor_visto / mayor_original
    resultado["multiplicadores"] = multiplicadores
    resultado["traza"] = traza
    return resultado

def sustitucion_regresiva(matriz_triangular, vector_nuevo):
    tamano = len(matriz_triangular)
    multdiv = 0         # contador de multiplicaciones y divisiones
    sumres = 0          # contador de sumas y restas

    incognitas = []                 # toda la x parte en 0
    for fila in range(tamano):
        incognitas.append(0.0)

    fila = tamano - 1               # se parte desde la ultima fila por eso el -1
    while fila >= 0:                # se sube hasta la primera
        suma = vector_nuevo[fila]

        for columna in range(fila + 1, tamano):     # resto las x que ya conozco
            suma = suma - matriz_triangular[fila][columna] * incognitas[columna]
            multdiv = multdiv + 1       # 1 multiplicacion
            sumres = sumres + 1         # 1 resta

        incognitas[fila] = suma / matriz_triangular[fila][fila]    # despejo la x de la fila
        multdiv = multdiv + 1           # 1 division

        fila = fila - 1             # subo una fila

    return incognitas, multdiv, sumres


def resolver_con_lu(multiplicadores, matriz_triangular, orden_filas, vector_otro):
    # resuelve con OTRO vector reutilizando L y U (sin volver a eliminar)
    tamano = len(matriz_triangular)
    flops = 0               # cuento lo que cuesta una pasada

    vector_ordenado = []        # pongo el vector en el mismo orden de filas
    for fila in range(tamano):
        vector_ordenado.append(vector_otro[orden_filas[fila]])

    for fila in range(tamano):      # sustitucion hacia adelante con L
        for columna in range(fila):
            vector_ordenado[fila] = vector_ordenado[fila] - multiplicadores[fila][columna] * vector_ordenado[columna]
            flops = flops + 2       # 1 multiplicacion + 1 resta

    incognitas, multdiv, sumres = sustitucion_regresiva(matriz_triangular, vector_ordenado)
    flops = flops + multdiv + sumres

    return incognitas, flops

def resolver(matriz, vector, con_pivoteo):
    # junta todo: elimina, despeja y suma el total de flops
    resultado = eliminacion_gaussiana(matriz, vector, con_pivoteo)
    incognitas, multdiv, sumres = sustitucion_regresiva(resultado["matriz_triangular"], resultado["vector_nuevo"])
    resultado["incognitas"] = incognitas
    resultado["multdiv_regresiva"] = multdiv
    resultado["sumres_regresiva"] = sumres

    total = multdiv + sumres
    for paso in range(len(matriz) - 1):
        total = total + resultado["multdiv_por_paso"][paso] + resultado["sumres_por_paso"][paso]
    resultado["total_flops"] = total
    return resultado

def escribir_traza(resultado, titulo, archivo):
    archivo.write("=" * 80 + "\n")
    archivo.write(titulo + "\n")
    archivo.write("=" * 80 + "\n")
    numero_paso = 1
    for foto_matriz, foto_vector, fila_original in resultado["traza"]:
        archivo.write("\n--- Paso " + str(numero_paso) + " (pivote = fila original " + str(fila_original) + ")\n")
        for fila in range(len(foto_matriz)):
            linea = ""
            for columna in range(len(foto_matriz)):
                linea = linea + f"{foto_matriz[fila][columna]: .16e} "
            linea = linea + "| " + f"{foto_vector[fila]: .16e}"
            archivo.write(linea + "\n")
        numero_paso = numero_paso + 1
    archivo.write("\n")
    

if __name__ == "__main__":
    matriz, vector = cargar_datos()

    sin_pivoteo = resolver(matriz, vector, False)
    con_pivoteo = resolver(matriz, vector, True)

    # guardo la traza de las dos corridas en resultados/
    carpeta_resultados = os.path.join(carpeta_proyecto, "resultados")
    os.makedirs(carpeta_resultados, exist_ok=True)
    archivo = open(os.path.join(carpeta_resultados, "traza_eliminacion.txt"), "w")
    escribir_traza(sin_pivoteo, "Corrida (i): sin pivoteo", archivo)
    escribir_traza(con_pivoteo, "Corrida (ii): pivoteo parcial escalado", archivo)
    archivo.close()

    # pruebo que L y U sirven para resolver de nuevo sin eliminar otra vez
    x_otra_vez, flops_pasada = resolver_con_lu(con_pivoteo["multiplicadores"], con_pivoteo["matriz_triangular"], con_pivoteo["orden_filas"], vector)