import os

carpeta_proyecto = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def cargar_datos():
    matriz = []
    archivo = open("datos/matriz_a.csv")
    for linea in archivo:
        numeritos = linea.strip().split(",")
        filita = []
        for numerito in numeritos:
            filita.append(float(numerito))
        matriz.append(filita)
    archivo.close()

    vector = []
    archivo = open("datos/vector_b.csv")
    for linea in archivo:
        vector.append(float(linea.strip()))
    archivo.close()

    return matriz, vector


if __name__ == "__main__":
    matriz, vector = cargar_datos()
    print("Filas de la matriz:", len(matriz))
    print("a(1,1) =", matriz[0][0])
    print("b(1) =", vector[0])