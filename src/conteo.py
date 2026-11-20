import numpy as np

def calcular_cotas_teoricas(n):
    """Calcula las cotas teóricas de referencia para una matriz densa nxn."""
    minimo_teorico_a = (2.0 / 3.0) * (n ** 3)  
    cota_referencia_guia = 1486                # Cota de referencia indicada en T1
    return minimo_teorico_a, cota_referencia_guia

def resumir_estrategia(resultado_dict):
    """Extrae y consolida los contadores desde el diccionario devuelto por resolver() en gauss.py."""
    md_elim = np.sum(resultado_dict["multdiv_por_paso"])
    sr_elim = np.sum(resultado_dict["sumres_por_paso"])
    total_elim = md_elim + sr_elim  #Suma de multiplicaciones/divisiones y sumas/restas en la fase de eliminación gaussiana

    md_reg = resultado_dict["multdiv_regresiva"]
    sr_reg = resultado_dict["sumres_regresiva"]
    total_reg = md_reg + sr_reg     #Suma de multiplicaciones/divisiones y sumas/restas en la fase de sustitución regresiva

    total_md = md_elim + md_reg
    total_sr = sr_elim + sr_reg
    total_flops = resultado_dict["total_flops"] #Guarda el total de operaciones de punto flotante (flops) calculadas en la función resolver().

    # Sobrecostos de estabilidad (pivoteo parcial escalado)
    comp_busqueda = np.sum(resultado_dict["comparaciones_por_paso"])
    comp_escalas = resultado_dict["comparaciones_escalas"]
    total_comparaciones = comp_busqueda + comp_escalas
    div_escala = np.sum(resultado_dict["divisiones_escala_por_paso"])
    total_intercambios = np.sum(resultado_dict["intercambios_por_paso"])

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

if __name__ == "__main__":
    n = 12
    minimo_teorico_a, cota_referencia_guia = calcular_cotas_teoricas(n)
    print(f"Para n={n}:")
    print(f"Minimo teorico A: {minimo_teorico_a}")
    print(f"Cota de referencia guia: {cota_referencia_guia}")