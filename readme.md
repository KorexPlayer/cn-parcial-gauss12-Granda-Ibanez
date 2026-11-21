# El pivote — Gauss 12×12
Resuelve un sistema lineal de 12 incógnitas con eliminación gaussiana, contando las operaciones y cuidando la mantisa.

## Autoría
- Catalina Ibañez
- Carlos Granda

Computación Numérica (INF-321), Ingeniería Civil Informática, UCM, 2026

## Requisitos
- Python (3.12.4 o superior)
- Instalar dependencias:
```
pip install -r requirements.txt
```

## Ejecución
Desarrollamos una carpeta raíz para el proyecto:
```
python src/main.py
```
Que entrega .csv y .txt correspondientes: 

| Archivo generado | Qué es lo qué contiene |
|---|---|
| `resultados/solucion.csv` | las 12 incógnitas con todos sus dígitos |
| `resultados/traza_eliminacion.txt` | la matriz después de cada uno de los 11 pasos |
| `resultados/conteo_flops.md` | tabla de operaciones por paso |
| `resultados/residuo.txt` | r = b − Ax y su norma infinito |

## Estructura
```
cn-parcial-gauss12-Granda-Ibanez/
├─ datos/        # matriz A y vector b sin truncar
├─ src/          # gauss, conteo, verificación y main
├─ resultados/   # todo lo que genera el ejecutable main.py
├─ informe.md    # análisis completo
└─ readme.md     # este mismo archivo
```

## Resultado resumido

| Incógnita | Valor |
|---|---:|
| x₁ | 0.7384915627 |
| x₂ | −1.2059374168 |
| x₃ | 0.0043871952 |
| x₄ | −0.0917264538 |
| x₅ | 2.8461093725 |
| x₆ | −0.0062815409 |
| x₇ | 1.4738205916 |
| x₈ | −0.0384927164 |
| x₉ | 0.9152648037 |
| x₁₀ | −2.6071394852 |
| x₁₁ | 0.0518736294 |
| x₁₂ | −1.0849371625 |

| Total de flops | ‖r‖∞ |
|---:|---:|
| 1354 | 2.60×10⁻¹⁰ |

Solución con pivoteo parcial escalado y una pasada de refinamiento iterativo.

## Informe
El análisis completo está en [informe.md](informe.md).

## Uso
Trabajo académico para la asignatura INF-321, Computación númerica UCM. Solo uso educativo.
