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
 
### Objetivo
Resolver un sistema lineal denso de 12 incógnitas con eliminación gaussiana, cumpliendo dos cosas que compiten entre sí:
 
- **Costo mínimo:** no pasarse del mínimo teórico de la eliminación densa.
- **Máxima exactitud:** elegir el pivote y el orden de las operaciones para perder la menor cantidad posible de bits de mantisa.
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
La escala de cada fila es su número más grande en valor absoluto, s(i) = max |a(i,j)|. Lo calcula la función `calcular_escalas` de `src/gauss.py`.
 
| Fila | Escala s(i) | Tipo de fila |
|---|---:|---|
| F1 | 9.27×10⁻⁵ | pequeña |
| F2 | 9.96×10⁴ | grande |
| F3 | 9.64×10⁻³ | pequeña |
| F4 | 8.38×10³ | grande |
| F5 | 8.05×10⁻⁵ | pequeña |
| F6 | 9.67×10⁵ | grande |
| F7 | 6.33×10⁻² | pequeña |
| F8 | 8.30×10⁴ | grande |
| F9 | 9.96×10⁻⁵ | pequeña |
| F10 | 8.87×10³ | grande |
| F11 | 9.74×10⁻³ | pequeña |
| F12 | 9.84×10⁴ | grande |
 
Las filas impares son "pequeñas" (entre 10⁻⁵ y 10⁻²) y las pares son "grandes" (entre 10³ y 10⁶). Entre F5 y F6 hay casi **diez órdenes de magnitud** de diferencia.
 
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
 
- **Un multiplicador por fila.** m(i,k) se calcula una sola vez y se usa en toda la fila y también en b. Si se recalculara por columna, el costo subiría a ≈ 1750 flops.
- **No tocamos lo que ya es cero.** El barrido parte en la columna k+1. El elemento a(i,k) no se calcula como a(i,k) − m·a(k,k), porque por cómo se construye m ese resultado es cero. Por eso lo ponemos en `0.0` directo. Así además quedan ceros **exactos** y no un "casi cero" por redondeo.
- **Sin inversa y una sola factorización.** Nunca calculamos A⁻¹. Guardamos L (los multiplicadores), U (la matriz triangular) y el orden de las filas, y con eso la función `resolver_con_lu` resuelve con **otro** vector b sin volver a eliminar. Eso lo usamos en el refinamiento iterativo (sección 8), y cada pasada cuesta solo **276 flops**.
### Por qué pivoteo *escalado* y no solo parcial
El pivoteo parcial normal elige el mayor |a(i,k)|, y en esta matriz siempre ganarían las filas grandes (10⁵) solo por estar escritas en una escala mayor. El escalado divide por s(i), o sea, compara cada candidato **con el tamaño de su propia fila**. Así una fila de 10⁻⁵ no pierde injustamente frente a una de 10⁵.
 
## 4. Conteo de operaciones
 
**Convención:** 1 flop = 1 operación aritmética elemental (+, −, ×, ÷). Las comparaciones, las divisiones |a(i,k)|/s(i) del pivoteo y los intercambios de filas **no** son flops: se cuentan aparte como **sobrecosto de estabilidad**.
 
Los números salen de contadores que están **dentro del código**, justo al lado de cada operación (`multdiv_por_paso`, `sumres_por_paso`, etc. en `src/gauss.py`). No salen de una fórmula escrita a mano. La tabla completa también queda en `resultados/conteo_flops.md`.
 
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
 
### Cómo se lee un paso
En el paso k quedan r = 12 − k filas bajo el pivote. Cada una cuesta:
- 1 división (el multiplicador),
- r multiplicaciones y r restas (columnas k+1 a 12),
- 1 multiplicación y 1 resta (el vector b).
Por ejemplo, en el paso 1 (r = 11): mult/div = 11 × (1 + 11 + 1) = **143** y sum/res = 11 × (11 + 1) = **132**. 
 
### Comparación con la cota
 
| Estrategia | Flops de la cota | Nuestros flops | Diferencia | ¿Dentro de ±10 %? |
|---|---:|---:|---:|---|
| Eliminación densa bien hecha | ≈ 1486 | 1354 | −8.9 % | Sí |
| + pivoteo parcial escalado | ≈ 1486 + 66 comparaciones | 1354 + sobrecosto aparte | −8.9 % | Sí |
 
La diferencia con 1486 es exactamente **132 flops** (66 multiplicaciones + 66 restas). La cota del enunciado incluye calcular a(i,k) − m·a(k,k) en la columna k. Nosotros no lo calculamos porque sabemos que da cero y lo asignamos directo: 1354 + 132 = 1486. Los flops son **los mismos en las dos corridas**, porque el pivoteo no agrega aritmética de eliminación.
 
### Sobrecosto de estabilidad (contado aparte)
 
| Corrida | Comparaciones | Divisiones de escala | Intercambios |
|---|---:|---:|---:|
| Sin pivoteo | 0 | 0 | 0 |
| Pivoteo parcial escalado | 198 (132 para las escalas + 66 de búsqueda) | 77 | 8 |
 
Las 66 comparaciones de búsqueda coinciden con las "66 comparaciones" que da el enunciado. Las 132 extra son de calcular las 12 escalas s(i) una sola vez al inicio.
 
## 5. Análisis de mantisa y pivoteo
 
Corrimos **el mismo código** dos veces, cambiando solo `con_pivoteo`. El error se mide contra la solución exacta del sistema, calculada con aritmética racional (`fractions`) a partir de los datos originales en texto.
 
| Corrida | Multiplicador máximo | Factor de crecimiento ρ | Error relativo ‖x̂ − x‖∞ / ‖x‖∞ | Peor error por componente | Bits perdidos (peor componente) | Dígitos correctos (mínimo) |
|---|---:|---:|---:|---:|---:|---:|
| (i) Sin pivoteo | 2.75×10¹³ | 2.64×10³ | 2.08×10⁻¹⁰ | 7.51×10⁻⁸ | ≈ 29 | 7.1 |
| (ii) Pivoteo parcial escalado | 2.47×10⁹ | 1.78 | 2.36×10⁻¹³ | 8.54×10⁻¹¹ | ≈ 20 | 10.1 |
 
- **ρ** = (número más grande que apareció durante la eliminación) / (número más grande de A original).
- **Bits perdidos** = 53 − bits correctos, con bits correctos = −log₂(peor error por componente).
### Dígitos correctos de cada incógnita
 
| Incógnita | Sin pivoteo | Con pivoteo escalado |
|---|---:|---:|
| x₁ | 10.2 | 13.1 |
| x₂ | 9.3 | 12.3 |
| x₃ | 7.1 | 10.1 |
| x₄ | 8.6 | 11.6 |
| x₅ | 10.0 | 12.9 |
| x₆ | 7.1 | 10.1 |
| x₇ | 9.6 | 12.5 |
| x₈ | 8.8 | 11.8 |
| x₉ | 9.8 | 12.7 |
| x₁₀ | 9.9 | 12.9 |
| x₁₁ | 8.9 | 11.9 |
| x₁₂ | 9.6 | 12.5 |
 
Las peores siempre son x₃ y x₆, que son las incógnitas más chicas (≈ 10⁻³): un mismo error absoluto pesa mucho más en ellas en términos relativos.
 
### ¿Por qué el pivote a(1,1) es inaceptable?
Con a(1,1) = 1.3407×10⁻⁸ como pivote, los multiplicadores del paso 1 se disparan: el más grande es m(6,1) ≈ 2.75×10¹³. En la fila 6 se calcula a(6,j) − m·a(1,j), y el término m·a(1,j) llega a ≈ 2.5×10⁹, mientras que la fila 6 original vale como mucho 9.7×10⁵. Al alinear exponentes, los bits más bajos de la fila 6 se salen de los 53 bits de mantisa y **se pierden en ese solo paso, sin forma de recuperarlos**. Lo mismo pasa en las otras filas: en el paso 1 se pierden entre 10 y 13 bits por fila.
 
Después, la matriz crece unas 2.6×10³ veces (≈ 11 bits más) y el error se sigue acumulando. Al final se pierden ≈ 29 bits en el peor componente y solo quedan **7 dígitos correctos** en x₃ y x₆. Lo peligroso es que el sistema "se resuelve igual", porque ningún pivote es cero, y el resultado se ve razonable aunque tenga la mitad de los dígitos malos.
 
### ¿Qué hace el pivoteo escalado?
En la columna 1, la fila 1 tiene |a(1,1)|/s(1) ≈ 1.45×10⁻⁴ (la peor de todas) y la fila 8 tiene razón 1.0 (su elemento de la columna 1 es el más grande de su fila). Por eso el primer pivote es **F8**. El orden completo de pivotes queda:
 
| Paso | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Fila original usada | F8 | F2 | F1 | F6 | F9 | F11 | F12 | F5 | F4 | F10 | F3 | F7 |
 
La fila 1 se usa recién en el paso 3, cuando ya no compite con filas de 10⁵ en su columna. La fila 7 (la casi dependiente) queda **al final**.
 
### ¿Y el multiplicador de 2.47×10⁹ con pivoteo?
Parece grande, pero no es un problema. El escalado a veces elige como pivote una fila pequeña (de 10⁻⁵) y la usa para eliminar en una fila grande (de 10⁵). Ese multiplicador es solo la razón entre las escalas de las dos filas, así que m·a(k,j) queda del mismo tamaño que la fila que lo recibe. Lo que importa es que la matriz casi no crece (**ρ = 1.78**), y por eso la pérdida de bits queda controlada.
 
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
3. **El residuo decide.** La cota se construye con el peor caso; lo que de verdad dice cuán buena es la solución es el residuo y la comparación con la solución exacta (sección 8).
### Casi dependencia de la fila 7
Ajustando F7 como combinación de F3 y F5 por mínimos cuadrados:
 
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
Con el residuo exacto resolvemos **A·d = r reutilizando L y U** (`resolver_con_lu`, sin volver a factorizar, **276 flops** por pasada), y corregimos x ← x + d.
 
| | Antes | Después de 1 pasada |
|---|---:|---:|
| Error relativo ‖x̂ − x‖∞/‖x‖∞ | 2.36×10⁻¹³ | 6.0×10⁻¹⁷ |
| Dígitos correctos (mínimo) | 10.1 | 16.1 |
| ‖r‖∞ | 2.29×10⁻¹⁰ | 2.60×10⁻¹⁰ |
 
El error bajó unas **4000 veces**, pero el residuo **no bajó**. Esto no significa que empeoró: ‖r‖∞ lo domina la fila 6, que tiene escala 10⁶. Ahí, el solo hecho de guardar x en `float64` ya produce un residuo de ≈ s(6)·‖x‖·ε ≈ 10⁻¹⁰. O sea, llegamos al **piso de la doble precisión** y el residuo no puede bajar más aunque x sea casi exacta. Por eso el residuo solo, sin mirar el error, puede engañar.
 
### Control con biblioteca
`numpy.linalg.solve` (pivoteo parcial sin escalar) da un error relativo de 2.55×10⁻¹³ (`residuo.txt`), parecido al nuestro (2.36×10⁻¹³). Lo usamos solo como control: la solución reportada es la de nuestra implementación.
 
## 9. Conclusiones
 
**¿Qué compromiso costo/exactitud elegimos y por qué?**
Pivoteo parcial escalado + una pasada de refinamiento iterativo. La eliminación cuesta **1354 flops** (−8.9 % de la cota): un multiplicador por fila, sin volver a la subdiagonal y sin inversa. La estabilidad cuesta 198 comparaciones, 77 divisiones de escala y 8 intercambios, y no agrega ningún flop a la eliminación. A cambio, el peor dígito pasa de **7 a 10**, y con una pasada de refinamiento (276 flops más) llegamos al límite del `float64` (error 6.0×10⁻¹⁷). La secuencia más barata (sin pivoteo) es justamente la que pierde ≈ 29 bits.
 
**¿Qué pasó con el pivote minúsculo?**
a(1,1) = 1.34×10⁻⁸ genera multiplicadores de 10¹³ y destruye entre 10 y 13 bits por fila en el primer paso. No se puede aceptar aunque no sea cero: el sistema se resuelve igual, pero con la mitad de los dígitos malos.
 
**¿Qué haríamos distinto con n = 1000?**
- La eliminación costaría ≈ 2n³/3 ≈ 6.7×10⁸ flops, y el pivoteo solo ≈ n² ≈ 10⁶ comparaciones (≈ 0.15 %). Pivotear seguiría siendo "gratis" en proporción.
- Usaríamos una librería optimizada (LU por bloques de LAPACK) en vez de bucles de Python, porque con listas sería lentísimo.
- Equilibraríamos las filas una vez al inicio en vez de dividir por s(i) en cada búsqueda.
- Estimaríamos κ con un método barato (tipo Hager) en vez de calcular todos los valores singulares.
- La traza completa de cada paso ya no sería práctica: guardaríamos solo ρ, los pivotes y el multiplicador máximo por paso.
-- El residuo en precisión extendida con `Decimal` sería lento para n = 1000; lo calcularíamos con doble-doble o con precisión mixta.
## 10. Reproducibilidad
 
| Cosa | Versión |
|---|---|
| Python | >=3.10 (`python --version`) |
| numpy | >=1.24.0 |
| Comando | `pip install -r requirements.txt` y luego `python src/main.py` |
 
`python src/main.py` genera todo `resultados/`: `solucion.csv`, `traza_eliminacion.txt`, `conteo_flops.md` y `residuo.txt`.
 
### Uso de IA
Usamos Claude (Anthropic) como apoyo para ___ (por ejemplo: organizar el trabajo, revisar el código y redactar partes del informe). Todo el código fue escrito, ejecutado y revisado por nosotros, y los números del informe salen de correr el programa.
 
### Fuentes
- Enunciado de la Prueba Parcial «El pivote», Computación Numérica (INF-321), UCM, noviembre de 2026.
- R. L. Burden, J. D. Faires. *Análisis numérico*. Cap. 6 (eliminación gaussiana y estrategias de pivoteo).
- Material de clases de Computación Numérica (sistemas lineales, pivoteo y condicionamiento).
 