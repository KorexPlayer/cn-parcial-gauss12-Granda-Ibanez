import os
import math
from decimal import Decimal, getcontext
from fractions import Fraction

# Configuramos 50 dígitos de precisión decimal para el residuo y el refinamiento
getcontext().prec = 50

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

def calcular_residuo_extendido(A_str, b_str, x):
    """
    Calcula r = b - A*x usando Decimal de 50 dígitos y retorna r y su norma infinito.
    """
    n = len(A_str)
    r_dec = []
    for i in range(n):
        suma_ax = Decimal(0)
        for j in range(n):
            # Convertimos float x[j] vía str() para capturar su valor exacto en punto flotante
            suma_ax += Decimal(A_str[i][j]) * Decimal(str(x[j]))
        ri = Decimal(b_str[i]) - suma_ax
        r_dec.append(ri)

    norma_inf = max(abs(val) for val in r_dec)
    return r_dec, norma_inf

def calcular_error_relativo_inf(x_aprox, x_exacta):
    """Calcula ||x_aprox - x_exacta||_inf / ||x_exacta||_inf en precisión extendida."""
    diffs = [abs(Decimal(str(x_aprox[i])) - (Decimal(x_exacta[i].numerator) / Decimal(x_exacta[i].denominator))) for i in range(len(x_aprox))]
    normas_x = [abs(Decimal(x_exacta[i].numerator) / Decimal(x_exacta[i].denominator)) for i in range(len(x_aprox))]
    return max(diffs) / max(normas_x)

def contar_digitos_significativos(x_val, x_frac):
    """Determina cuántos dígitos significativos son correctos comparando con la solución exacta."""
    x_exacta_dec = Decimal(x_frac.numerator) / Decimal(x_frac.denominator)
    x_aprox_dec = Decimal(str(x_val))

    if x_exacta_dec == 0:
        return 16 if x_aprox_dec == 0 else 0

    err_rel = abs(x_aprox_dec - x_exacta_dec) / abs(x_exacta_dec)
    if err_rel == 0:
        return 16
    # Cantidad de dígitos significativos correctos: -log10(error_relativo)
    digitos = math.floor(-float(err_rel.ln() / Decimal(10).ln()))
    return max(0, min(16, digitos))


