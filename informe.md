# El pivote — Eliminación gaussiana con mínimo costo y máxima mantisa
Prueba parcial Unidad nro 3 computación numérica.

<!-- ESQUELETO: todo lo que está entre ___ o entre <!-- --> hay que reemplazarlo o borrarlo. -->

## 1. Identificación

| Campo | Detalle |
|---|---|
| Integrantes | Catalina Ibañez · Carlos ___ |
| Asignatura | Computación Numérica (INF-321) — Ingeniería Civil Informática, UCM |
| Docentes | Juan Manuel Cabezas Contardo · Pedro Pesante |
| Fecha | ___ de noviembre de 2026 |
| Repositorio | https://github.com/KorexPlayer/cn-parcial-gauss12-Granda-Ibanez |

### Objetivo
Resolver un sistema lineal denso de 12 incógnitas con eliminación gaussiana, cumpliendo dos cosas que compiten entre sí:

- **Costo mínimo:** no pasarse del mínimo teórico de flops (≈ 2n³/3).
- **Máxima exactitud:** elegir el pivote y el orden de las operaciones para perder la menor cantidad de bits de mantisa.

Nota sobre la convención de conteo: 1 flop = 1 operación aritmética elemental (+, −, ×, ÷). Las comparaciones e intercambios del pivoteo se cuentan aparte. ___

## 2. Planteamiento

### Ambiente controlado

| Indicador | Valor |
|---|---|
| Tamaño del sistema | n = 12 |
| Elementos no nulos | 144 de 144 (densa) |
| Cifras por elemento de A | 9 |
| Dígitos de b | 20 |
| Rango de magnitudes en A | 10⁻⁸ … 10⁵ |
| det(A) | ≈ 1.3646×10⁹ (no singular) |
| κ₂(A) | ≈ 3.08×10¹³ |
| Precisión usada | doble IEEE 754 (53 bits, ε ≈ 2.22×10⁻¹⁶) |

### Escalas de fila

| Fila | Orden de magnitud | Comentario |
|---|---|---|
| F1 | 10⁻⁵ (a(1,1) = 10⁻⁸) | ___ |
| F2 | 10⁴ | ___ |
| F3 | 10⁻³ | ___ |
| F4 | 10³ | ___ |
| F5 | 10⁻⁵ | ___ |
| F6 | 10⁵ | ___ |
| F7 | 10⁻² | casi combinación de F3 y F5 |
| F8 | 10⁴ | ___ |
| F9 | 10⁻⁵ | ___ |
| F10 | 10³ | ___ |
| F11 | 10⁻³ | ___ |
| F12 | 10⁴ | ___ |

### Qué lo hace difícil
<!-- Explicar con tus palabras los 3 problemas: -->
1. **Escalas de fila incompatibles:** ___
2. **Pivote inicial minúsculo:** ___
3. **Casi dependencia lineal oculta (fila 7):** ___

## 3. Estrategia de eliminación

| Paso | Qué se hace | Fórmula |
|---|---|---|
| 1 | Elegir pivote (escalado) | max \|a(i,k)\| / s(i), con s(i) = max_j \|a(i,j)\| |
| 2 | Calcular el multiplicador una sola vez | m(i,k) = a(i,k) / a(k,k) |
| 3 | Actualizar la fila desde la columna k+1 | a(i,j) = a(i,j) − m·a(k,j) |
| 4 | Actualizar b con el mismo multiplicador | b(i) = b(i) − m·b(k) |
| 5 | Sustitución regresiva | x(i) = (b(i) − Σ a(i,j)x(j)) / a(i,i) |

### Por qué no nos pasamos del mínimo
<!-- multiplicador una vez por fila, barrido desde k+1, sin inversa, L y U reutilizadas -->
___

## 4. Conteo de operaciones

Datos de `resultados/conteo_flops.md` (salen del contador del código, no de una fórmula).

| Paso k | Mult/Div | Sum/Res | Flops del paso | Comparaciones | Intercambios |
|---:|---:|---:|---:|---:|---:|
| 1 | ___ | ___ | ___ | ___ | ___ |
| 2 | ___ | ___ | ___ | ___ | ___ |
| 3 | ___ | ___ | ___ | ___ | ___ |
| 4 | ___ | ___ | ___ | ___ | ___ |
| 5 | ___ | ___ | ___ | ___ | ___ |
| 6 | ___ | ___ | ___ | ___ | ___ |
| 7 | ___ | ___ | ___ | ___ | ___ |
| 8 | ___ | ___ | ___ | ___ | ___ |
| 9 | ___ | ___ | ___ | ___ | ___ |
| 10 | ___ | ___ | ___ | ___ | ___ |
| 11 | ___ | ___ | ___ | ___ | ___ |
| Sustitución regresiva | ___ | ___ | ___ | — | — |
| **Total** | ___ | ___ | ___ | ___ | ___ |

### Comparación con la cota

| Estrategia | Flops cota | Flops nuestros | Diferencia | ¿Dentro de ±10 %? |
|---|---:|---:|---:|---|
| Eliminación densa bien hecha | ≈ 1486 | ___ | ___ % | ___ |
| + pivoteo parcial | ≈ 1486 + 66 comparaciones | ___ | ___ | ___ |

___

## 5. Análisis de mantisa y pivoteo

### Resultados de las dos corridas

| Corrida | Multiplicador máximo | Factor de crecimiento | Error relativo | Bits perdidos | Dígitos correctos |
|---|---:|---:|---:|---:|---:|
| (i) Sin pivoteo | ___ | ___ | ___ | ___ | ___ |
| (ii) Pivoteo parcial escalado | ___ | ___ | ___ | ___ | ___ |

### ¿Por qué el pivote a(1,1) es inaceptable?
___

### ¿Por qué difieren las dos corridas (en bits de mantisa)?
___

### Orden de pivotes con escalado

| Paso | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Fila original usada | ___ | ___ | ___ | ___ | ___ | ___ | ___ | ___ | ___ | ___ | ___ |

## 6. Condicionamiento

| Dato | Valor |
|---|---|
| κ(A) estimado | ___ |
| Cota κ·ε | ___ |
| Dígitos que garantiza la cota | ___ |
| Dígitos obtenidos de verdad | ___ |

### ¿Por qué la cota es pesimista aquí?
___

### Casi dependencia de la fila 7
F7 ≈ ___ · F3 + ___ · F5 (perturbación de orden ___)

## 7. Resultados

| Incógnita | Valor reportado | Dígitos correctos afirmados |
|---|---|---:|
| x₁ | ___ | ___ |
| x₂ | ___ | ___ |
| x₃ | ___ | ___ |
| x₄ | ___ | ___ |
| x₅ | ___ | ___ |
| x₆ | ___ | ___ |
| x₇ | ___ | ___ |
| x₈ | ___ | ___ |
| x₉ | ___ | ___ |
| x₁₀ | ___ | ___ |
| x₁₁ | ___ | ___ |
| x₁₂ | ___ | ___ |

## 8. Verificación

| Prueba | Resultado |
|---|---|
| ‖r‖∞ = ‖b − Ax‖∞ (precisión extendida) | ___ |
| ‖r‖∞ después de 1 refinamiento | ___ |
| Diferencia máx. con `numpy.linalg.solve` (solo control) | ___ |

___

## 9. Conclusiones

**¿Qué compromiso costo/exactitud elegimos y por qué?** ___

**¿Qué pasó con el pivote minúsculo?** ___

**¿Qué haríamos distinto con n = 1000?** ___

## 10. Reproducibilidad

| Cosa | Versión |
|---|---|
| Python | ___ |
| numpy | ___ |
| mpmath | ___ |
| Comando | `python src/main.py` |

### Uso de IA
Usamos ___ para ___. Podemos explicar cada línea del código.

### Fuentes
- ___
