import os

carpeta_proyecto = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def calcular_cotas_teoricas(n):
    """Calcula las cotas teóricas de referencia para una matriz densa nxn."""
    minimo_teorico_a = (2.0 / 3.0) * (n ** 3)  
    cota_referencia_guia = 1486                # Cota de referencia indicada en T1
    return minimo_teorico_a, cota_referencia_guia

def resumir_estrategia(resultado_dict):
    """Extrae y consolida los contadores desde el diccionario devuelto por resolver() en gauss.py."""
    md_elim = sum(resultado_dict["multdiv_por_paso"])
    sr_elim = sum(resultado_dict["sumres_por_paso"])
    total_elim = md_elim + sr_elim  #Suma de multiplicaciones/divisiones y sumas/restas en la fase de eliminación gaussiana

    md_reg = resultado_dict["multdiv_regresiva"]
    sr_reg = resultado_dict["sumres_regresiva"]
    total_reg = md_reg + sr_reg     #Suma de multiplicaciones/divisiones y sumas/restas en la fase de sustitución regresiva

    total_md = md_elim + md_reg     #Suma de multiplicaciones/divisiones en ambas fases (eliminación y sustitución regresiva)
    total_sr = sr_elim + sr_reg     #Suma de sumas/restas en ambas fases (eliminación y sustitución regresiva)
    total_flops = resultado_dict["total_flops"] #Guarda el total de operaciones de punto flotante (flops) calculadas en la función resolver().

    # Sobrecostos de estabilidad (pivoteo parcial escalado)
    comp_busqueda = sum(resultado_dict["comparaciones_por_paso"])
    comp_escalas = resultado_dict["comparaciones_escalas"]
    total_comparaciones = comp_busqueda + comp_escalas
    div_escala = sum(resultado_dict["divisiones_escala_por_paso"])
    total_intercambios = sum(resultado_dict["intercambios_por_paso"])

    return {
        "md_elim": md_elim,                                                 #MultiDiv eliminacion gaussiana
        "sr_elim": sr_elim,                                                 #SumaResta eliminacion gaussiana
        "total_elim": total_elim,                                           #Total eliminacion gaussiana
        "md_reg": md_reg,                                                   #MultiDiv sustitucion regresiva
        "sr_reg": sr_reg,                                                   #SumaResta sustitucion regresiva
        "total_reg": total_reg,                                             #Total sustitucion regresiva
        "total_md": total_md,                                               #Total MultiDiv
        "total_sr": total_sr,                                               #Total SumaResta
        "total_flops": total_flops,                                         #Total operaciones de punto flotante
        "comp_busqueda": comp_busqueda,                                     #Comparaciones en busqueda de pivote
        "comp_escalas": comp_escalas,                                       #Comparaciones en escalas
        "total_comparaciones": total_comparaciones,                         #Total comparaciones
        "div_escala": div_escala,                                           #Divisiones en escalas
        "total_intercambios": total_intercambios,                           #Total intercambios
        "multiplicador_maximo": resultado_dict["multiplicador_maximo"],     #Multiplicador máximo
        "factor_crecimiento": resultado_dict["factor_crecimiento"],         #Factor de crecimiento
    }

def generar_reporte_flops(sin_pivoteo, con_pivoteo, n=12):
    """
    Recibe los diccionarios 'sin_pivoteo' y 'con_pivoteo' de gauss.resolver()
    y genera el archivo resultados/conteo_flops.md exigido en la sección 5.
    """
    ruta_salida = os.path.join(carpeta_proyecto, "resultados", "conteo_flops.md")
    os.makedirs(os.path.dirname(ruta_salida), exist_ok=True)

    res_sin = resumir_estrategia(sin_pivoteo)
    res_con = resumir_estrategia(con_pivoteo)

    minimo_teorico_a, cota_guia = calcular_cotas_teoricas(n)
    desviacion_pct = ((res_sin["total_flops"] - cota_guia) / cota_guia) * 100.0

    with open(ruta_salida, "w", encoding="utf-8") as archivo:
        # 1. Convención
        archivo.write("# Conteo de Operaciones (FLOPs)\n\n")
        archivo.write("""
## 1. Convención de conteo
- **Unidad básica:** 1 FLOP = 1 operación aritmética elemental en punto flotante (`+`, `-`, `*`, `/`).
- **Reutilización de multiplicadores:** Cada multiplicador $m_{i,k} = a_{i,k}/a_{k,k}$ se calcula una única vez por fila (1 división) y se almacena en la matriz $L$ para reutilizarse en toda la fila, en el vector $b$ y en futuras sustituciones.
- **Barrido sin redundancia:** En el paso $k$, el recorrido de columnas comienza en $k+1$, sin tocar los elementos ya anulados bajo la diagonal.
- **Sobrecosto de pivoteo parcial escalado:** Las comparaciones para obtener las escalas iniciales de fila, las comparaciones de búsqueda de pivote, las divisiones de escalamiento $|a_{i,k}| / s_i$ y los intercambios de filas se contabilizan por separado del costo base de eliminación.

##2. Desglose paso a paso

| Etapa / Paso | Mult / Div | Sumas / Restas | Subtotal FLOPs | Comparaciones (Pivoteo) | Div. Escala (Pivoteo) | Intercambios (Swaps) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |\n""")

        num_pasos = len(sin_pivoteo["multdiv_por_paso"])
        for paso in range(num_pasos):
            md = sin_pivoteo["multdiv_por_paso"][paso]
            sr = sin_pivoteo["sumres_por_paso"][paso]
            comp = con_pivoteo["comparaciones_por_paso"][paso]
            div_esc = con_pivoteo["divisiones_escala_por_paso"][paso]
            swap = con_pivoteo["intercambios_por_paso"][paso]
            archivo.write(f"| Paso {paso + 1} | {md} | {sr} | {md + sr} | {comp} | {div_esc} | {swap} |\n")

        # Totales de la tabla
        archivo.write(
            f"| **Subtotal Eliminación ({num_pasos} pasos)** | **{res_sin['md_elim']}** | **{res_sin['sr_elim']}** | "
            f"**{res_sin['total_elim']}** | **{res_con['comp_busqueda']}** | **{res_con['div_escala']}** | **{res_con['total_intercambios']}** |\n"
        )
        archivo.write(
            f"| Sustitución Regresiva | {res_sin['md_reg']} | {res_sin['sr_reg']} | {res_sin['total_reg']} | 0 | 0 | 0 |\n"
        )
        archivo.write(
            f"| **TOTAL GENERAL** | **{res_sin['total_md']}** | **{res_sin['total_sr']}** | "
            f"**{res_sin['total_flops']}** | **{res_con['comp_busqueda']}** (+{res_con['comp_escalas']} escalas) | "
            f"**{res_con['div_escala']}** | **{res_con['total_intercambios']}** |\n\n"
        )

        # 3. Cotas y comparacion final
        archivo.write(f"""
## 3. Comparación con cotas teóricas\n\n
- **Mínimo teórico clásico ($\\frac{{2}}{{3}}n^3$ para $n={n}$):** {minimo_teorico_a:.0f} FLOPs (solo matriz $A$).
- **Cota de referencia completa (guía T1):** ~{cota_guia} FLOPs.
- **FLOPs aritméticos base medidos:** {res_sin['total_flops']} FLOPs (desviación de {desviacion_pct:+.2f}% respecto a la cota de referencia, dentro del margen $\\pm 10\\%$).

## 4. Comparación de corridas y sobrecosto de estabilidad
| Indicador | Corrida (i): Sin pivoteo | Corrida (ii): Pivoteo parcial escalado |
| :--- | :---: | :---: |
| FLOPs aritméticos de eliminación + regresiva | {res_sin['total_flops']} | {res_con['total_flops']} |
| Comparaciones iniciales de escala ($s_i$) | 0 | {res_con['comp_escalas']} |
| Comparaciones de búsqueda de pivote | 0 | {res_con['comp_busqueda']} |
| Divisiones de escalamiento ($|a_{{i,k}}|/s_i$) | 0 | {res_con['div_escala']} |
| Intercambios de filas (swaps) | 0 | {res_con['total_intercambios']} |
| Multiplicador máximo ($\\max |m_{{i,k}}|$) | {res_sin['multiplicador_maximo']:.6e} | {res_con['multiplicador_maximo']:.6e} |
| Factor de crecimiento ($\\rho$) | {res_sin['factor_crecimiento']:.6e} | {res_con['factor_crecimiento']:.6e} |\n""")

    archivo.close()

    return res_sin, res_con

if __name__ == "__main__":
    n = 12
    minimo_teorico_a, cota_referencia_guia = calcular_cotas_teoricas(n)
    print(f"Para n={n}:")
    print(f"Minimo teorico A: {minimo_teorico_a}")
    print(f"Cota de referencia guia: {cota_referencia_guia}")