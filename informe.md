# El pivote — Eliminación gaussiana con mínimo costo y máxima mantisa
Prueba parcial Unidad nro 3 computación numérica.
 
## 1. Identificación
 
| Campo | Detalle |
|---|---|
| Integrantes | Catalina Ibañez y Carlos Granda |
| Asignatura | Computación Numérica (INF-321) — Ingeniería Civil Informática, UCM |
| Docente | Juan Manuel Cabezas Contardo |
| Sección | 2 |
| Fecha | Talca, 20 de noviembre de 2026 |
| Repositorio | https://github.com/KorexPlayer/cn-parcial-gauss12-Granda-Ibanez |
 
## 2. Planteamiento
 
En este trabajo resolvemos un sistema lineal **Ax = b** de 12 ecuaciones y 12 incógnitas usando eliminación gaussiana. La matriz A es **densa**: sus 144 elementos son diferentes de cero, por lo que no existen ceros que nos ayuden a reducir operaciones. Cada elemento de A tiene 9 cifras significativas y el vector b viene con 20 dígitos. Todos se guardaron completos en `datos/matriz_a.csv` y `datos/vector_b.csv`, sin truncar.
 
Trabajamos en doble precisión IEEE 754, que tiene 53 bits de mantisa (ε ≈ 2.22×10⁻¹⁶).
 
### Ambiente controlado
 
| Indicador | Valor |
|---|---|
| Tamaño del sistema | n = 12 |
| Elementos no nulos | 144 de 144 |
| Cifras por elemento de A | 9 |
| Dígitos de b | 20 |
| Rango de magnitudes en A | 10⁻⁸ … 10⁶ |
| det(A) | ≈ 1.3646×10⁹ (no singular) |
| κ₂(A) | ≈ 3.08×10¹³ |
 
### Escalas de fila
La escala de cada fila es s(i) = max |a(i,j)| (`calcular_escalas` en `src/gauss.py`). Las filas impares son "pequeñas" (s(i) entre 10⁻⁵ y 10⁻²) y las pares son "grandes" (entre 10³ y 10⁶). Entre F5 (8.05×10⁻⁵) y F6 (9.67×10⁵) hay casi **diez órdenes de magnitud**.
 
### Qué hace difícil este sistema
 
1. **Escalas de fila incompatibles.** El computador debe igualar los exponentes cuando a una fila pequeña se le resta un múltiplo de una fila grande. En ese procedimiento, los bits más bajos del número pequeño se escapan de los 53 bits de la mantisa y se extravían para siempre.
2. **Pivote inicial minúsculo.** El primer pivote es a(1,1) = 1.3407×10⁻⁸. No solo es pequeño comparado con su columna (a(2,1) ≈ 1.19×10⁴), sino también dentro de su propia fila: a(1,1)/s(1) ≈ 1.45×10⁻⁴. Si lo usamos como pivote, los multiplicadores se vuelven gigantes.
3. **Casi dependencia lineal oculta.** La fila 7 está muy cerca de ser una combinación de las filas 3 y 5. Por eso el número de condición es tan alto (κ₂ ≈ 3.08×10¹³) aunque la matriz no sea singular (sección 6).
Por todo esto, el costo de este problema no está solo en los flops, sino también en los **bits de mantisa** que se pueden perder. Para ver ese efecto, comparamos dos corridas del mismo código: una **sin pivoteo** y otra con **pivoteo parcial escalado**.

 
## 3. Estrategia de eliminación
 
Todo está en `src/gauss.py`. Usamos **el mismo código** para las dos corridas; solo cambia el parámetro `con_pivoteo` (`False` o `True`) de la función `resolver`.
 
### Pasos del algoritmo
 
| Paso | Qué se hace | Fórmula | Función |
|---|---|---|---|
| 1 | Calcular la escala de cada fila (una sola vez, solo con pivoteo) | s(i) = max_j \|a(i,j)\| | `calcular_escalas` |
| 2 | Elegir el pivote | fila con mayor \|a(i,k)\| / s(i) | `elegir_fila_pivote` |
| 3 | Intercambiar filas si hace falta (en A, b, L, escalas y orden de filas) | — | `eliminacion_gaussiana` |
| 4 | Calcular el multiplicador **una sola vez por fila** y guardarlo en L | m(i,k) = a(i,k) / a(k,k) | `eliminacion_gaussiana` |
| 5 | Actualizar la fila desde la columna k+1 | a(i,j) = a(i,j) − m·a(k,j) | `eliminacion_gaussiana` |
| 6 | Actualizar b con el mismo multiplicador | b(i) = b(i) − m·b(k) | `eliminacion_gaussiana` |
| 7 | Sustitución regresiva | x(i) = (b(i) − Σ a(i,j)·x(j)) / a(i,i) | `sustitucion_regresiva` |
 
Se hacen **11 pasos** (k = 1, …, 11) y quedan los **66 elementos bajo la diagonal en cero**. Después de cada paso guardamos una "foto" de la matriz y del vector, y todas quedan en `resultados/traza_eliminacion.txt`.
 
### Por qué no nos pasamos del mínimo
 
- **Un multiplicador por fila**, reutilizado en toda la fila y en b (recalcularlo por columna costaría ≈ 1750 flops).
- **No tocamos lo que ya es cero:** el barrido parte en k+1 y a(i,k) se asigna 0.0 directo, así quedan ceros exactos.
- **Sin inversa y una sola factorización:** guardamos L, U y el orden de filas; `resolver_con_lu` los reutiliza en el refinamiento (276 flops por pasada).
 
## 4. Conteo de operaciones
 
**Convención:** 1 flop = 1 operación aritmética elemental (+, −, ×, ÷). Las comparaciones, las divisiones |a(i,k)|/s(i) del pivoteo y los intercambios de filas **no** son flops: se cuentan aparte como **sobrecosto de estabilidad**.
 
Los números salen de contadores dentro de `src/gauss.py`, al lado de cada operación (no de una fórmula). La tabla completa queda en `resultados/conteo_flops.md`.
 
| Paso k | Filas bajo el pivote | Mult/Div | Sum/Res | Flops del paso |
|---:|---:|---:|---:|---:|
| 1 | 11 | 143 | 132 | 275 |
| 2 | 10 | 120 | 110 | 230 |
| 3 | 9 | 99 | 90 | 189 |
| 4 | 8 | 80 | 72 | 152 |
| 5 | 7 | 63 | 56 | 119 |
| 6 | 6 | 48 | 42 | 90 |
| 7 | 5 | 35 | 30 | 65 |
| 8 | 4 | 24 | 20 | 44 |
| 9 | 3 | 15 | 12 | 27 |
| 10 | 2 | 8 | 6 | 14 |
| 11 | 1 | 3 | 2 | 5 |
| **Subtotal eliminación (incluye b)** | | **638** | **572** | **1210** |
| Sustitución regresiva | | 78 | 66 | 144 |
| **Total** | | **716** | **638** | **1354** |
 
Por ejemplo, en el paso 1 hay 11 filas bajo el pivote: 11 × (1 + 11 + 1) = 143 mult/div y 11 × (11 + 1) = 132 sum/res. **Total: 1354 flops, −8.9 % de la cota ≈ 1486 (dentro del ±10 %).**
 
La diferencia con 1486 es exactamente **132 flops** (66 multiplicaciones + 66 restas). La cota del enunciado incluye calcular a(i,k) − m·a(k,k) en la columna k. Nosotros no lo calculamos porque sabemos que da cero y lo asignamos directo: 1354 + 132 = 1486. Los flops son **los mismos en las dos corridas**, porque el pivoteo no agrega aritmética de eliminación.
 
### Sobrecosto de estabilidad (contado aparte)
 
| Corrida | Comparaciones | Divisiones de escala | Intercambios |
|---|---:|---:|---:|
| Sin pivoteo | 0 | 0 | 0 |
| Pivoteo parcial escalado | 198 (132 para las escalas + 66 de búsqueda) | 77 | 8 |
 
 
## 5. Análisis de mantisa y pivoteo
 
Corrimos **el mismo código** dos veces, cambiando solo `con_pivoteo`. El error se mide contra la solución exacta del sistema, calculada con aritmética racional (`fractions`) a partir de los datos originales en texto.
 
| Corrida | Multiplicador máximo | Factor de crecimiento ρ | Error relativo ‖x̂ − x‖∞ / ‖x‖∞ | Peor error por componente | Bits perdidos (peor componente) | Dígitos correctos (mínimo) |
|---|---:|---:|---:|---:|---:|---:|
| (i) Sin pivoteo | 2.75×10¹³ | 2.64×10³ | 2.08×10⁻¹⁰ | 7.51×10⁻⁸ | ≈ 29 | 7.1 |
| (ii) Pivoteo parcial escalado | 2.47×10⁹ | 1.78 | 2.36×10⁻¹³ | 8.54×10⁻¹¹ | ≈ 20 | 10.1 |
 
- **ρ** = (número más grande que apareció durante la eliminación) / (número más grande de A original).
- **Bits perdidos** = 53 − bits correctos, con bits correctos = −log₂(peor error por componente).

Las peores siempre son x₃ y x₆, que son las incógnitas más chicas (≈ 10⁻³): un mismo error absoluto pesa mucho más en ellas en términos relativos.
 
### ¿Por qué el pivote a(1,1) es inaceptable?
Con a(1,1) = 1.3407×10⁻⁸ como pivote, los multiplicadores del paso 1 se disparan: el más grande es m(6,1) ≈ 2.75×10¹³. En la fila 6 se calcula a(6,j) − m·a(1,j), y el término m·a(1,j) llega a ≈ 2.5×10⁹, mientras que la fila 6 original vale como mucho 9.7×10⁵. Al alinear exponentes, los bits más bajos de la fila 6 se salen de los 53 bits de mantisa y **se pierden en ese solo paso, sin forma de recuperarlos**. Lo mismo pasa en las otras filas: en el paso 1 se pierden entre 10 y 13 bits por fila.
 
Después, la matriz crece unas 2.6×10³ veces (≈ 11 bits más) y el error se sigue acumulando. Al final se pierden ≈ 29 bits en el peor componente y solo quedan **7 dígitos correctos** en x₃ y x₆. Lo peligroso es que el sistema "se resuelve igual", porque ningún pivote es cero, y el resultado se ve razonable aunque tenga la mitad de los dígitos malos.
 
### ¿Qué hace el pivoteo escalado?
Compara |a(i,k)|/s(i): en la columna 1, F1 tiene 1.45×10⁻⁴ (la peor) y F8 tiene 1.0, así que el primer pivote es **F8**. Orden de pivotes: F8, F2, F1, F6, F9, F11, F12, F5, F4, F10, F3, F7. La fila 7 (la casi dependiente) queda al final. Aunque el multiplicador máximo es 2.47×10⁹, es solo la razón entre escalas de fila: la matriz casi no crece (ρ = 1.78).
 
### Resumen de la comparación
El pivoteo escalado no cambia los flops (1354 en las dos), pero reduce el factor de crecimiento de 2641 a 1.78 y el error relativo unas **880 veces** (de 2.08×10⁻¹⁰ a 2.36×10⁻¹³). Se pasa de 7 a 10 dígitos correctos en la peor incógnita: **≈ 3 dígitos ≈ 10 bits** recuperados, a cambio de 198 comparaciones, 77 divisiones de escala y 8 intercambios.
 
## 6. Condicionamiento
 
| Dato | Valor |
|---|---|
| κ₂(A) (con `numpy.linalg.cond`, solo como control) | 3.08×10¹³ |
| Cota κ₂·ε | 6.85×10⁻³ |
| Dígitos que garantiza la cota | 2.16 |
| Dígitos obtenidos con pivoteo (sin refinar) | entre 10.1 y 13.1 |
 
La cota clásica dice que ‖δx‖/‖x‖ ≲ κ(A)·ε ≈ 6.85×10⁻³, o sea que **no garantiza ni 2 dígitos**. Pero en la práctica obtuvimos **10 dígitos** como mínimo con pivoteo escalado.
 
### ¿Por qué la cota es pesimista aquí?
1. **Es el peor caso.** La cota supone que los errores de redondeo se acomodan justo en la dirección que más se amplifica, y eso casi nunca pasa.
2. **Gran parte de κ viene de las escalas, no de la matriz "en sí".** Si dividimos cada fila por su escala s(i) (equilibrar filas), el número de condición baja de 3.08×10¹³ a **≈ 8.1×10³**. Con ese κ, la cota da ≈ 8.1×10³ · 2.22×10⁻¹⁶ ≈ 1.8×10⁻¹² → **≈ 11–12 dígitos**, que calza con lo que obtuvimos. El pivoteo escalado trabaja justamente sobre esa parte.
### Casi dependencia de la fila 7
Ajustando F7 como combinación de F3 y F5 por mínimos cuadrados (cálculo hecho aparte con numpy; no forma parte de `main.py`):
 
**F7 ≈ 4.177 · F3 + 583.25 · F5**, con un residuo relativo de **2.3×10⁻³**.
 
O sea, F7 está a una perturbación chiquita de ser combinación lineal de F3 y F5. Eso hace que la matriz esté "casi" singular. El pivoteo escalado lo deja en evidencia solo: usa F7 como **último pivote**, y ese pivote final vale U(12,12) ≈ 1.57×10⁻⁴, apenas 2.5×10⁻³ de la escala de su fila (s(7) = 6.33×10⁻²), del mismo orden que el residuo del ajuste.
 
## 7. Resultados
 
Valores de `resultados/solucion.csv`, corrida con pivoteo parcial escalado. La columna "sin refinar" se escribe solo con los dígitos que afirmamos correctos (comparados contra la solución exacta con `fractions`). La columna "refinado" es tras una pasada de refinamiento iterativo, con error relativo de 6.0×10⁻¹⁷ (sección 8).

| Incógnita | Valor sin refinar | Dígitos afirmados | Valor refinado |
|---|---|---:|---|
| x₁ | 0.7384915626999 | 13 | 0.7384915627 |
| x₂ | −1.2059374168 | 12 | −1.2059374168 |
| x₃ | 0.0043871952 | 10 | 0.0043871952 |
| x₄ | −0.0917264538 | 11 | −0.0917264538 |
| x₅ | 2.8461093725 | 12 | 2.8461093725 |
| x₆ | −0.006281540899 | 10 | −0.0062815409 |
| x₇ | 1.4738205916 | 12 | 1.4738205916 |
| x₈ | −0.0384927164 | 11 | −0.0384927164 |
| x₉ | 0.9152648037 | 12 | 0.9152648037 |
| x₁₀ | −2.6071394852 | 12 | −2.6071394852 |
| x₁₁ | 0.0518736294 | 11 | 0.0518736294 |
| x₁₂ | −1.0849371625 | 12 | −1.0849371625 |

Todas superan los 8 dígitos mínimos. Las peores son x₃ y x₆ (10 dígitos), las incógnitas más chicas.
 
## 8. Verificación
 
### Residuo en precisión extendida
Calculamos r = b − Ax en **precisión extendida** (`Decimal` con 50 dígitos), usando A y b desde el texto original, para que el residuo no se contamine con el redondeo del `float`. La solución exacta de referencia se calcula con `fractions`.
 
| Corrida | ‖r‖∞ |
|---|---:|
| Sin pivoteo | 1.77×10⁻⁷ |
| Con pivoteo escalado (antes de refinar) | 2.29×10⁻¹⁰ |
| Con pivoteo escalado (después de 1 refinamiento) | 2.60×10⁻¹⁰ |
 
### Refinamiento iterativo
Con ese residuo resolvemos **A·d = r reutilizando L y U** (`resolver_con_lu`, sin volver a factorizar, **276 flops** por pasada), y corregimos x ← x + d.
 
| | Antes | Después de 1 pasada |
|---|---:|---:|
| Error relativo ‖x̂ − x‖∞/‖x‖∞ | 2.36×10⁻¹³ | 6.0×10⁻¹⁷ |
| Dígitos correctos (mínimo) | 10.1 | 16.1 |
 
El error bajó unas **4000 veces**, pero el residuo no: ‖r‖∞ lo domina la fila 6 (escala 10⁶), donde solo guardar x en `float64` ya da ≈ 10⁻¹⁰. Llegamos al piso de la doble precisión.
 
### Control con biblioteca
`numpy.linalg.solve` (pivoteo parcial sin escalar) da un error relativo de 2.55×10⁻¹³ (`residuo.txt`), parecido al nuestro (2.36×10⁻¹³). Lo usamos solo como control: la solución reportada es la de nuestra implementación.
 
## 9. Conclusiones
 
**¿Qué compromiso costo/exactitud elegimos y por qué?**
Pivoteo parcial escalado + una pasada de refinamiento iterativo. La eliminación cuesta **1354 flops** (−8.9 % de la cota): un multiplicador por fila, sin volver a la subdiagonal y sin inversa. La estabilidad cuesta 198 comparaciones, 77 divisiones de escala y 8 intercambios, y no agrega ningún flop a la eliminación. A cambio, el peor dígito pasa de **7 a 10**, y con una pasada de refinamiento (276 flops más) llegamos al límite del `float64` (error 6.0×10⁻¹⁷). La secuencia más barata (sin pivoteo) es justamente la que pierde ≈ 29 bits.
 
 
**¿Qué haríamos distinto con n = 1000?**
- Costaría ≈ 6.7×10⁸ flops y el pivoteo solo ≈ 10⁶ comparaciones (0.15 %), así que pivotear seguiría siendo barato.
- Usaríamos LU por bloques de una librería optimizada (LAPACK) en vez de bucles de Python, y equilibraríamos las filas una vez al inicio.
- Guardaríamos solo ρ, los pivotes y el multiplicador máximo por paso, no la traza completa.

## 10. Reproducibilidad
 
| Cosa | Versión |
|---|---|
| Python | 3.12.4 |
| numpy | 2.3.1 |
| Comando | `pip install -r requirements.txt` y luego `python src/main.py` |
 
`python src/main.py` genera todo `resultados/`: `solucion.csv`, `traza_eliminacion.txt`, `conteo_flops.md` y `residuo.txt`.
 
### Uso de IA
Usamos Claude (Anthropic) como apoyo para organizar el trabajo, guiar la escritura del código paso a paso y redactar el informe. Todo el código fue escrito, ejecutado y revisado por nosotros, y los números del informe salen de correr el programa.
 
### Fuentes
- Enunciado de la Prueba Parcial «El pivote», Computación Numérica (INF-321), UCM, noviembre de 2026.
- R. L. Burden, J. D. Faires. *Análisis numérico*. Cap. 6 (eliminación gaussiana y estrategias de pivoteo).
- Material de clases de Computación Numérica (sistemas lineales, pivoteo y condicionamiento).
 