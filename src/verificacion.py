import os
from decimal import Decimal, getcontext
from fractions import Fraction


carpeta_proyecto = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def cargar_datos_exactos():
    """Lee matriz_a.csv y vector_b.csv directamente como texto para convertirlos sin pérdida."""
    ruta_a = os.path.join(carpeta_proyecto, "datos", "matriz_a.csv")
    ruta_b = os.path.join(carpeta_proyecto, "datos", "vector_b.csv")

    A_str = []
    with open(ruta_a, "r", encoding="utf-8") as f:
        for linea in f:
            A_str.append([val.strip() for val in linea.strip().split(",")])

    b_str = []
    with open(ruta_b, "r", encoding="utf-8") as f:
        for linea in f:
            b_str.append(linea.strip())

    return A_str, b_str


def resolver_exacto_fractions(A_str, b_str):
    """
    Resuelve Ax = b en aritmética racional exacta (Fraction).
    Sirve como referencia verdadera para medir el error relativo (T2)
    y los dígitos significativos correctos (T3).
    """
    n = len(A_str)
    A = [[Fraction(val) for val in fila] for fila in A_str]
    b = [Fraction(val) for val in b_str]

    for k in range(n - 1):
        # Pivoteo simple para evitar división por cero en racionales
        if A[k][k] == 0:
            for i in range(k + 1, n):
                if A[i][k] != 0:
                    A[k], A[i] = A[i], A[k]
                    b[k], b[i] = b[i], b[k]
                    break

        for i in range(k + 1, n):
            factor = A[i][k] / A[k][k]
            for j in range(k + 1, n):
                A[i][j] -= factor * A[k][j]
            b[i] -= factor * b[k]

    x_exacta = [Fraction(0)] * n
    for i in range(n - 1, -1, -1):
        suma = b[i]
        for j in range(i + 1, n):
            suma -= A[i][j] * x_exacta[j]
        x_exacta[i] = suma / A[i][i]

    return x_exacta