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

def eliminacion_gaussiana(matriz, vector):
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

    multdiv_por_paso = []       # contador de multiplicaciones y divisiones
    sumres_por_paso = []        # contador de sumas y restas
    for paso in range(tamano - 1):
        multdiv_por_paso.append(0)
        sumres_por_paso.append(0)

    for fila_pivote in range(tamano - 1):       # cada vueltita usa una fila como pivote
        pivote = matriz_copia[fila_pivote][fila_pivote]

        for fila_de_abajo in range(fila_pivote + 1, tamano):

            numero_a_eliminar = matriz_copia[fila_de_abajo][fila_pivote]
            multiplicador = numero_a_eliminar / pivote
            multdiv_por_paso[fila_pivote] = multdiv_por_paso[fila_pivote] + 1      # conte 1 division

            matriz_copia[fila_de_abajo][fila_pivote] = 0.0      #este numero ya se queda en cero, lo pongo directo

            for columna in range(fila_pivote + 1, tamano):      #esto actualiza el resto de la fila
                numero_de_arriba = matriz_copia[fila_pivote][columna]
                matriz_copia[fila_de_abajo][columna] = matriz_copia[fila_de_abajo][columna] - multiplicador * numero_de_arriba
                multdiv_por_paso[fila_pivote] = multdiv_por_paso[fila_pivote] + 1   # 1 multiplicacion
                sumres_por_paso[fila_pivote] = sumres_por_paso[fila_pivote] + 1     # 1 resta

            vector_copia[fila_de_abajo] = vector_copia[fila_de_abajo] - multiplicador * vector_copia[fila_pivote]       #al vector le pasa exactamente lo mismo, con el mismo multiplicador
            multdiv_por_paso[fila_pivote] = multdiv_por_paso[fila_pivote] + 1       # 1 multiplicacion
            sumres_por_paso[fila_pivote] = sumres_por_paso[fila_pivote] + 1         # 1 resta

    return matriz_copia, vector_copia, multdiv_por_paso, sumres_por_paso

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

if __name__ == "__main__":
    matriz, vector = cargar_datos()
    matriz_triangular, vector_nuevo = eliminacion_gaussiana(matriz, vector)


    cantidad_de_ceros = 0       # cuento los ceros que hay
    for fila in range(12):
        for columna in range(fila):
            if matriz_triangular[fila][columna] == 0.0:
                cantidad_de_ceros = cantidad_de_ceros + 1

    print("Ceros bajo la diagonal:", cantidad_de_ceros)

    incognitas = sustitucion_regresiva(matriz_triangular, vector_nuevo)

    print("Solucion sin pivoteo:")
    for fila in range(12):
        print("x" + str(fila + 1) + " =", incognitas[fila])