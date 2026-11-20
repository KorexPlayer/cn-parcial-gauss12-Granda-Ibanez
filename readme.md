# El pivote — Gauss 12×12
Resuelve un sistema lineal denso de 12 incógnitas con eliminación gaussiana, contando las operaciones y cuidando la mantisa.

## Autoría
- Catalina Ibañez
- Carlos ___

Computación Numérica (INF-321) · Ingeniería Civil Informática · UCM · 2026

## Requisitos
- Python ___ (3.10 o superior)
- Instalar dependencias:
```
pip install -r requirements.txt
```

## Ejecución
Desde la carpeta raíz del proyecto:
```
python src/main.py
```

| Archivo generado | Qué tiene |
|---|---|
| `resultados/solucion.csv` | las 12 incógnitas con todos sus dígitos |
| `resultados/traza_eliminacion.txt` | la matriz después de cada uno de los 11 pasos |
| `resultados/conteo_flops.md` | tabla de operaciones por paso |
| `resultados/residuo.txt` | r = b − Ax y su norma infinito |

## Estructura
```
cn-parcial-gauss12-Granda-Ibanez/
├─ datos/        # matriz A y vector b sin truncar
├─ src/          # código: gauss, conteo, verificación y main
├─ resultados/   # todo lo que genera main.py
├─ informe.md    # análisis completo
└─ readme.md     # este archivo
```

## Resultado resumido

| Incógnita | Valor (10 dígitos) |
|---|---|
| x₁ | ___ |
| x₂ | ___ |
| x₃ | ___ |
| x₄ | ___ |
| x₅ | ___ |
| x₆ | ___ |
| x₇ | ___ |
| x₈ | ___ |
| x₉ | ___ |
| x₁₀ | ___ |
| x₁₁ | ___ |
| x₁₂ | ___ |

| Total de flops | ‖r‖∞ |
|---|---|
| ___ | ___ |

## Informe
El análisis completo está en [informe.md](informe.md).

## Uso
Trabajo académico para la asignatura INF-321, UCM. Solo uso educativo.
