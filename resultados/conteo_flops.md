# Conteo de Operaciones (FLOPs)


## 1. Convención de conteo
- **Unidad básica:** 1 FLOP = 1 operación aritmética elemental en punto flotante (`+`, `-`, `*`, `/`).
- **Reutilización de multiplicadores:** Cada multiplicador $m_{i,k} = a_{i,k}/a_{k,k}$ se calcula una única vez por fila (1 división) y se almacena en la matriz $L$ para reutilizarse en toda la fila, en el vector $b$ y en futuras sustituciones.
- **Barrido sin redundancia:** En el paso $k$, el recorrido de columnas comienza en $k+1$, sin tocar los elementos ya anulados bajo la diagonal.
- **Sobrecosto de pivoteo parcial escalado:** Las comparaciones para obtener las escalas iniciales de fila, las comparaciones de búsqueda de pivote, las divisiones de escalamiento $|a_{i,k}| / s_i$ y los intercambios de filas se contabilizan por separado del costo base de eliminación.

## 2. Desglose paso a paso

| Etapa / Paso | Mult / Div | Sumas / Restas | Subtotal FLOPs | Comparaciones (Pivoteo) | Div. Escala (Pivoteo) | Intercambios (Swaps) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| Paso 1 | 143 | 132 | 275 | 11 | 12 | 1 |
| Paso 2 | 120 | 110 | 230 | 10 | 11 | 0 |
| Paso 3 | 99 | 90 | 189 | 9 | 10 | 1 |
| Paso 4 | 80 | 72 | 152 | 8 | 9 | 1 |
| Paso 5 | 63 | 56 | 119 | 7 | 8 | 1 |
| Paso 6 | 48 | 42 | 90 | 6 | 7 | 1 |
| Paso 7 | 35 | 30 | 65 | 5 | 6 | 1 |
| Paso 8 | 24 | 20 | 44 | 4 | 5 | 1 |
| Paso 9 | 15 | 12 | 27 | 3 | 4 | 1 |
| Paso 10 | 8 | 6 | 14 | 2 | 3 | 0 |
| Paso 11 | 3 | 2 | 5 | 1 | 2 | 0 |
| **Subtotal Eliminación (11 pasos)** | **638** | **572** | **1210** | **66** | **77** | **8** |
| Sustitución Regresiva | 78 | 66 | 144 | 0 | 0 | 0 |
| **TOTAL GENERAL** | **716** | **638** | **1354** | **66** (+132 escalas) | **77** | **8** |


## 3. Comparación con cotas teóricas


- **Mínimo teórico clásico ($\frac{2}{3}n^3$ para $n=12$):** 1152 FLOPs (solo matriz $A$).
- **Cota de referencia completa (guía T1):** ~1486 FLOPs.
- **FLOPs aritméticos base medidos:** 1354 FLOPs (desviación de -8.88% respecto a la cota de referencia, dentro del margen $\pm 10\%$).

## 4. Comparación de corridas y sobrecosto de estabilidad

| Indicador | Corrida (i): Sin pivoteo | Corrida (ii): Pivoteo parcial escalado |
| :--- | :---: | :---: |
| FLOPs aritméticos de eliminación + regresiva | 1354 | 1354 |
| Comparaciones iniciales de escala ($s_i$) | 0 | 132 |
| Comparaciones de búsqueda de pivote | 0 | 66 |
| Divisiones de escalamiento ($|a_{i,k}|/s_i$) | 0 | 77 |
| Intercambios de filas (swaps) | 0 | 8 |
| Multiplicador máximo ($\max |m_{i,k}|$) | 2.752593e+13 | 2.465558e+09 |
| Factor de crecimiento ($\rho$) | 2.641224e+03 | 1.783202e+00 |
