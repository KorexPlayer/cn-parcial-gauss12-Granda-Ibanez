import os
import math
import csv
import numpy as np
import gauss as gs
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
    Sirve como referencia verdadera para medir el error relativo
    y los dígitos significativos correctos.
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
            suma_ax += Decimal(A_str[i][j]) * Decimal(float(x[j]))
        ri = Decimal(b_str[i]) - suma_ax
        r_dec.append(ri)

    norma_inf = max(abs(val) for val in r_dec)
    return r_dec, norma_inf

def calcular_error_relativo_inf(x_aprox, x_exacta):
    """Calcula ||x_aprox - x_exacta||_inf / ||x_exacta||_inf en precisión extendida."""
    diffs = [abs(Decimal(float(x_aprox[i])) - (Decimal(x_exacta[i].numerator) / Decimal(x_exacta[i].denominator))) for i in range(len(x_aprox))]
    normas_x = [abs(Decimal(x_exacta[i].numerator) / Decimal(x_exacta[i].denominator)) for i in range(len(x_aprox))]
    return max(diffs) / max(normas_x)

def contar_digitos_significativos(x_val, x_frac):
    """Determina cuántos dígitos significativos son correctos comparando con la solución exacta."""
    x_exacta_dec = Decimal(x_frac.numerator) / Decimal(x_frac.denominator)
    x_aprox_dec = Decimal(float(x_val))

    if x_exacta_dec == 0:
        return 16 if x_aprox_dec == 0 else 0

    err_rel = abs(x_aprox_dec - x_exacta_dec) / abs(x_exacta_dec)
    if err_rel == 0:
        return 16
    # Cantidad de dígitos significativos correctos: -log10(error_relativo)
    digitos = math.floor(-float(err_rel.ln() / Decimal(10).ln()))
    return max(0, min(16, digitos))

def ejecutar_verificacion(sin_pivoteo, con_pivoteo):
    """
    Ejecuta toda la tarea T4 y C5:
    1. Residuo en precisión extendida (sin y con pivoteo).
    2. Refinamiento iterativo reutilizando L y U (resolver_con_lu).
    3. Control con numpy.linalg.solve y condicionamiento kappa_2(A).
    4. Genera resultados/residuo.txt y resultados/solucion.csv.
    """
    A_str, b_str = cargar_datos_exactos()
    n = len(A_str)

    # Solución exacta en racionales para evaluación rigurosa de errores
    x_exacta = resolver_exacto_fractions(A_str, b_str)

    x_sin = sin_pivoteo["incognitas"]
    x_con = con_pivoteo["incognitas"]

    # 1. Residuos en precisión extendida
    r_sin, norm_r_sin = calcular_residuo_extendido(A_str, b_str, x_sin)
    r_con, norm_r_con = calcular_residuo_extendido(A_str, b_str, x_con)

    err_rel_sin = calcular_error_relativo_inf(x_sin, x_exacta)
    err_rel_con = calcular_error_relativo_inf(x_con, x_exacta)

    # 2. Refinamiento iterativo (1 pasada reutilizando L y U de la corrida con pivoteo)
    r_con_float = [float(val) for val in r_con]
    resultado_lu = gs.resolver_con_lu(
        con_pivoteo["multiplicadores"],
        con_pivoteo["matriz_triangular"],
        con_pivoteo["orden_filas"],
        r_con_float,
    )
    d = resultado_lu[0] if isinstance(resultado_lu, tuple) else resultado_lu
    x_refinada = [x_con[i] + d[i] for i in range(n)]
    r_ref, norm_r_ref = calcular_residuo_extendido(A_str, b_str, x_refinada)
    err_rel_ref = calcular_error_relativo_inf(x_refinada, x_exacta)

    # 3. Control con biblioteca (NumPy) y 4. Condicionamiento kappa(A)
    A_np = np.array([[float(val) for val in fila] for fila in A_str], dtype=np.float64)
    b_np = np.array([float(val) for val in b_str], dtype=np.float64)

    x_numpy = np.linalg.solve(A_np, b_np)
    r_np, norm_r_np = calcular_residuo_extendido(A_str, b_str, x_numpy)
    err_rel_np = calcular_error_relativo_inf(x_numpy, x_exacta)

    kappa_2 = np.linalg.cond(A_np, 2)
    eps_ieee = np.finfo(np.float64).eps  # ~2.220446e-16
    cota_error = kappa_2 * eps_ieee
    digitos_garantizados_cota = max(0.0, -math.log10(cota_error))

    # Guardar resultados/residuo.txt
    carpeta_resultados = os.path.join(carpeta_proyecto, "resultados")
    os.makedirs(carpeta_resultados, exist_ok=True)
    ruta_residuo = os.path.join(carpeta_resultados, "residuo.txt")

    with open(ruta_residuo, "w", encoding="utf-8") as f:
        f.write("=" * 80 + "\n")
        f.write("VERIFICACIÓN EN PRECISIÓN EXTENDIDA, REFINAMIENTO Y CONDICIONAMIENTO\n")
        f.write("=" * 80 + "\n\n")

        f.write("1. CONDICIONAMIENTO Y COTA DE DÍGITOS RECUPERABLES (C5)\n")
        f.write(f"   - Número de condición kappa_2(A) : {kappa_2:.6e}\n")
        f.write(f"   - Épsilon de máquina (IEEE 754)  : {eps_ieee:.6e}\n")
        f.write(f"   - Cota clásica kappa(A) * eps    : {cota_error:.6e}\n")
        f.write(f"   - Dígitos garantizados por cota  : {digitos_garantizados_cota:.2f} dígitos (cota pesimista)\n\n")

        f.write("2. COMPARACIÓN DE NORMAS DE RESIDUO ||r||_inf Y ERROR RELATIVO\n")
        f.write(f"   - Sin pivoteo          -> ||r||_inf = {float(norm_r_sin):.6e} | Error rel. = {float(err_rel_sin):.6e}\n")
        f.write(f"   - Con pivoteo escalado -> ||r||_inf = {float(norm_r_con):.6e} | Error rel. = {float(err_rel_con):.6e}\n")
        f.write(f"   - Tras 1 refinamiento  -> ||r||_inf = {float(norm_r_ref):.6e} | Error rel. = {float(err_rel_ref):.6e}\n")
        f.write(f"   - Control NumPy solve  -> ||r||_inf = {float(norm_r_np):.6e} | Error rel. = {float(err_rel_np):.6e}\n\n")

        f.write("3. VECTOR RESIDUO r = b - Ax (CON PIVOTEO PARCIAL ESCALADO, 50 DÍGITOS)\n")
        for i in range(n):
            f.write(f"   r[{i + 1:2d}] = {r_con[i]:+.15e}\n")
        f.write(f"\n   Norma infinito ||r||_inf = {norm_r_con:.15e}\n\n")

        f.write("4. VECTOR RESIDUO TRAS REFINAMIENTO ITERATIVO (1 PASADA)\n")
        for i in range(n):
            f.write(f"   r_ref[{i + 1:2d}] = {r_ref[i]:+.15e}\n")
        f.write(f"\n   Norma infinito ||r_ref||_inf = {norm_r_ref:.15e}\n")

    # Guardar resultados/solucion.csv con los dígitos significativos defendibles (T3 / C4)
    ruta_solucion = os.path.join(carpeta_resultados, "solucion.csv")
    with open(ruta_solucion, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["incognita", "valor_defendible", "digitos_correctos_afirmados", "valor_refinado"])
        for i in range(n):
            dig_con = contar_digitos_significativos(x_con[i], x_exacta[i])
            dig_ref = contar_digitos_significativos(x_refinada[i], x_exacta[i])
            # Formateamos la incógnita exactamente a los dígitos significativos defendibles
            val_formateado = f"{x_con[i]:.{dig_con}g}"
            val_ref_formateado = f"{x_refinada[i]:.{dig_ref}g}"
            writer.writerow([f"x{i + 1}", val_formateado, dig_con, val_ref_formateado])

    return {
        "kappa_2": kappa_2,
        "cota_error": cota_error,
        "norm_r_sin": float(norm_r_sin),
        "norm_r_con": float(norm_r_con),
        "norm_r_ref": float(norm_r_ref),
        "err_rel_sin": float(err_rel_sin),
        "err_rel_con": float(err_rel_con),
        "err_rel_ref": float(err_rel_ref),
        "x_refinada": x_refinada,
    }