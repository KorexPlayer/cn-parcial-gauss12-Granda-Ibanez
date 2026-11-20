
def calcular_cotas_teoricas(n):
    """Calcula las cotas teóricas de referencia para una matriz densa nxn."""
    minimo_teorico_a = (2.0 / 3.0) * (n ** 3)  
    cota_referencia_guia = 1486                # Cota de referencia indicada en T1
    return minimo_teorico_a, cota_referencia_guia

if __name__ == "__main__":
    n = 12
    minimo_teorico_a, cota_referencia_guia = calcular_cotas_teoricas(n)
    print(f"Para n={n}:")
    print(f"Minimo teorico A: {minimo_teorico_a}")
    print(f"Cota de referencia guia: {cota_referencia_guia}")