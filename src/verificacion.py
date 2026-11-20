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
