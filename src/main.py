import gauss as gs
import conteo as ct
import os
carpeta_proyecto = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    n = 12  # Tamaño del sistema
    A, b = gs.cargar_datos()

    resultado_sin_pivoteo = gs.resolver(A, b, con_pivoteo=False)

    resultado_con_pivoteo = gs.resolver(A, b, con_pivoteo=True)

    carpeta_resultados = os.path.join(carpeta_proyecto, "resultados")
    os.makedirs(carpeta_resultados, exist_ok=True)
    archivo = open(os.path.join(carpeta_resultados, "traza_eliminacion.txt"), "w")
    gs.escribir_traza(resultado_sin_pivoteo, "Corrida (i): sin pivoteo", archivo)
    gs.escribir_traza(resultado_con_pivoteo, "Corrida (ii): pivoteo parcial escalado", archivo)
    archivo.close()

    # Generar reporte de FLOPs
    ct.generar_reporte_flops(resultado_sin_pivoteo, resultado_con_pivoteo, n)

if __name__ == "__main__":
    main()