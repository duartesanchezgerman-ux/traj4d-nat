# Auditoría de los cálculos del artículo A1 (3 nov 2026; repetida con 70 días el 5 nov)

> **5 nov 2026, versión final con 70 días:** 61 comprobaciones, 0 fallos (`resultados/resultados_auditoria_a1.txt`). La sección J ya no lleva cifras fijas: reconstruye cada frase del manuscrito con los datos y la busca literal en `main.tex`. Las cifras finales están en el apartado 6. Los apartados 1–5 describen la auditoría del 3 nov con 10 días; los métodos no han cambiado.

Este documento explica, paso a paso, de dónde sale cada número del artículo y cómo se ha comprobado. Cada comprobación se repite con un comando, sin fiarse de la versión anterior.

**Comprobación automática:** `PYTHONPATH=src python scripts/auditoria_a1.py`
Resultado: `resultados/resultados_auditoria_a1.txt`, con **0 fallos en 49 comprobaciones**. El script tiene secciones A–J; abajo se cita cada una por su letra.

La auditoría encontró **tres defectos**. Los tres están corregidos y el análisis se ha rehecho entero:
1. **Hora local en vez de UTC** en el volcado de OpenSky (ver `12-correccion-hora-adsc.md`).
2. **Informes fechados una hora tarde en los extremos del vuelo.** La primera corrección de la vuelta de hora solo los detectaba en mitad del vuelo; ahora también al principio y al final: 1009 informes en vez de 971.
3. **Se podía perder el último informe del tramo.** Al quitar informes repetidos a menos de 5 km, la ruta volada podía acabar hasta 4 km antes que la ortodrómica y el óptimo con los que se compara. Afectaba a 13 de 376 tramos y rebajaba un poco el exceso. Corregido en `traj4d.ruta.puntos_distintos`, que conserva siempre el primer y el último informe.

**Cifras finales del A1 (10 días):**
- Exceso medio de la ruta volada frente al óptimo: **0,30 % (IC95 0,19–0,44)**.
- Mediana: 0,13 %. Unos 26 kg por tramo.
- Parte estructural (la rejilla de puntos de paso): 0,063 %, el **21 %** del exceso.
- Error del viento de ERA5 frente al medido a bordo: **3,11 m/s**.

---

## 1. Datos de entrada

| Dato | De dónde | Cómo se comprobó |
|---|---|---|
| Informes ADS-C (621 183 con posición) | Volcado de OpenSky en Zenodo | El recuento coincide con el que da la descripción del conjunto de datos. El decodificador solo usa la biblioteca estándar. El Mac y Claude obtienen tablas idénticas (sha256 en `docs/ADSC_UTC_SHA256.txt`). |
| Hora de cada informe | Campo Timestamp, en hora local de Europa central | 444 informes cruzados con ADS-B (adsb.lol, UTC). Tras convertir a UTC, la diferencia mediana es de 1 s (sección D y `resultados_f2_hora_adsc.txt`). |
| Tipo de avión | Base de aeronaves de OpenSky | El 98,8 % de los informes del Atlántico Norte tiene tipo. |
| Viento y temperatura | ERA5 por horas, 150–350 hPa, 0,25° | Comparación con lo que miden los aviones (sección RQ2 del artículo). |
| Modelo de consumo | Poll–Schumann (pycontrails 0.65.0) | Validado frente a EUROCONTROL/EEA (sección RQ3; la CERT de la OACI se retiró el 5 nov). |

## 2. Cadena de cálculo y comprobación de cada paso

### 2.1 Atmósfera estándar (sección A)
- **Qué hace:** el nivel de vuelo (FL) se pasa a presión con la atmósfera estándar ISA. Es lo que mide el altímetro a 1013,25 hPa.
- **Comprobación:** las presiones de FL300, 350, 370, 390 y 410 coinciden con la tabla ISA hasta la centésima de hPa. Por ejemplo, FL390 da 196,77 hPa.

### 2.2 Viento medido por el avión (sección B)
- **Qué hace:** el ADS-C da la dirección DE la que sopla el viento. Las componentes son u = −V·sin(dir) hacia el este y v = −V·cos(dir) hacia el norte.
- **Comprobación:** con este convenio, el error frente a ERA5 es de 3,1 m/s. Con el contrario sería de 73,9 m/s, imposible con un viento medio de 32 m/s. El convenio es el correcto.

### 2.3 Interpolación de ERA5 (sección C)
- **Qué hace:** interpolación lineal en tiempo, logaritmo de la presión, latitud y longitud.
- **Comprobación:** 300 informes al azar, recalculados con otra herramienta (xarray), dan lo mismo hasta 10⁻¹³.
- **Límite honesto:** xarray usa por debajo la misma biblioteca numérica (scipy). Lo que comprueba esta prueba es el manejo de las coordenadas: unidades de tiempo, orden de las latitudes, signo de las longitudes y paso de FL a presión. Ahí es donde suelen estar los errores. La comprobación física de que ERA5 está bien leído es la 2.2: ERA5 y los aviones casan a 3 m/s.

### 2.4 Hora y posiciones (sección D)
- **Velocidad sobre el suelo medida a bordo (etiqueta 14) frente a distancia entre informes / tiempo entre informes:** diferencia mediana de +0,13 kt en 111 tramos. Hora y posición son coherentes entre sí.
- **Duración del tramo calculada por el modelo** (ruta volada, Mach, viento ERA5) **frente a la observada:** diferencia mediana de −0,26 % con Mach medido y −0,22 % con el Mach del tipo. Cuadra.
- **Hallazgo:** el código mide las distancias sobre una esfera (R = 6371 km). Sobre el elipsoide real (WGS84) los tramos son un 0,32 % más largos.
  - Esto rebaja en esa proporción el combustible y el tiempo absolutos.
  - En el exceso relativo, que compara rutas a latitudes muy parecidas, el efecto medio es de +0,0002 puntos porcentuales: despreciable.
  - Se declara como limitación (ix) en el artículo. Para el A2, donde importan los tiempos absolutos, conviene pasar a WGS84.

### 2.5 Consumo (sección E)
- **Comprobación de orden de magnitud:** consumo horario medio en crucero a 0,8·MTOW, por tipo. Por ejemplo, B789 5411 kg/h, A359 5755 kg/h, B77W 8182 kg/h. Todos caen en el rango típico público de crucero.
- **Límite honesto:** esos rangos son aproximados y los he puesto yo a partir de valores públicos generales. Sirven para detectar un error grosero (un factor 2 o un cambio de unidades), no para validar el modelo. La validación de verdad es la de la sección RQ3 del artículo, frente a EUROCONTROL/EEA.

### 2.6 Integrador (sección F)
- **Qué hace el código:** Runge–Kutta de orden 4 con pasos de 20 km a lo largo de la ruta. El viento se toma en la posición y la hora de cada paso, y la masa baja con el consumo.
- **Comprobación:** un integrador distinto, escrito aparte (Euler, pasos de 1 km), en 25 tramos. La diferencia media es de +0,008 % y la máxima de 0,23 %. Esa máxima sale de que una versión interpola en ortodrómica y la otra en línea recta en latitud-longitud.

### 2.7 Optimizador (sección G)
- **Qué hace:** programación dinámica sobre una rejilla de 41 rutas paralelas a la ortodrómica, separadas 15 km (±300 km). La ruta elegida se vuelve a evaluar con el integrador completo.
- **Comprobaciones:**
  - El óptimo nunca es peor que la ortodrómica (peor caso −0,003 %).
  - En 61 de 376 tramos la ruta volada sale ligeramente mejor que el "óptimo". Es la resolución de la rejilla: mediana −0,04 %, peor caso −0,12 %, todos dentro del suelo de ruido declarado de ±0,2 %.
  - Antes de corregir el defecto 3, el peor caso era −0,32 %. Ese caso salía precisamente del último informe perdido.

### 2.8 Descomposición estructura × elección (sección H)
- **Qué hace:** separa el exceso en dos partes.
  - Lo que cuesta volar por la rejilla de puntos de paso (latitudes enteras o medias en los meridianos de 10°).
  - Lo que cuesta haber elegido unos puntos de esa rejilla y no los mejores.
- **Identidad:** (1 + exceso total) = (1 + estructural) × (1 + elección). Se cumple en todos los tramos con error de 4·10⁻¹⁶.
- **Coherencia:** el óptimo restringido nunca es mejor que el libre más allá del ruido (peor caso −0,10 %).

### 2.9 Estadística (sección I)
- **Qué hace:** bootstrap por bloques de días. Se remuestrean días, no tramos, porque los tramos de un mismo día comparten el tiempo meteorológico y no son independientes.
- **Comprobación:** otra implementación con otra semilla y 4 veces más remuestras da 0,30 % [0,19; 0,44]. Es lo mismo que en el artículo.

### 2.10 Cifras del manuscrito (sección J)
24 cifras clave del texto y de las tablas en inglés se recalculan desde las tablas por tramo y se comparan con el texto literal. Todas coinciden. El texto en español usa las mismas cifras. Esta comprobación no se automatizó para el español; se revisó a mano.

## 3. Qué cambió en el A1 con todas las correcciones

| Cifra | Versión del 2 nov | Final (3 nov) |
|---|---|---|
| Error del viento de ERA5 (RMS) | 4,28 m/s | 3,11 m/s |
| Infravaloración con viento > 60 m/s | −1,53 m/s | −0,86 m/s |
| Tramos evaluados | 379 | 376 |
| Exceso medio | 0,30 % (0,19–0,44) | 0,30 % (0,19–0,44) |
| Parte estructural | 0,056 % (19 %) | 0,063 % (21 %) |
| Pruebas automáticas | 23 | 28 |

**Conclusión:** el resultado principal se mantiene y la validación de ERA5 sale bastante mejor. Los tres defectos estaban en la preparación de los datos, no en el modelo ni en el optimizador.

## 4. Cómo puedes verificarlo tú (Germán)
1. Abre `resultados/resultados_auditoria_a1.txt`: cada línea dice [OK] o [FALLO] y qué se comprueba.
2. Para repetirlo en el Mac, cuando estén los 70 días: el script 2 rehace todo. Después lanza `PYTHONPATH=src python scripts/auditoria_a1.py`. Con 70 días las cifras de la sección J cambiarán y habrá que actualizar el texto; el resto de secciones debe seguir en [OK].
3. Tres comprobaciones sencillas a mano:
   - **Unidades.** Un tramo medio dura 81 min y quema unos 8570 kg (unos 6350 kg/h). 26 kg de exceso entre 8570 kg es un 0,30 %, igual que el exceso medio del artículo.
   - **Hora.** En `datos/procesado/hora_adsc_vs_adsb.csv`, la columna `desfase_si_local_s` vale −7200 en verano y −3600 en invierno.
   - **Tabla 1, consistencia interna.** La media ponderada de los cuatro grupos de estación y sentido debe dar la media global: (130 × 0,434 + 17 × 0,160 + 83 × 0,245 + 146 × 0,228) / 376 = 0,300 %. En la tabla aparecen redondeados (0,43; 0,16; 0,24; 0,23).

## 5. Limitaciones que siguen abiertas (declaradas en el artículo)
- Solo 10 días: la extensión a 70 está en marcha y es la que debe ir al artículo.
- Mach medido en el propio vuelo solo en el 22 % de los tramos.
- Nivel de vuelo constante: el A1 no evalúa la optimización vertical.
- Masa desconocida, aunque el resultado no depende de ella: 0,204 / 0,203 / 0,202 % a 70 / 80 / 90 % de la MTOW.
- Distancias sobre la esfera (efecto despreciable en el exceso relativo).
- ERA5 asimila observaciones de aviones: su error es una cota inferior del error de una previsión.
- Algunos informes fechados una hora tarde pueden seguir sin detectarse: 1 de 444 en la comprobación con ADS-B, al final de un vuelo con los informes anteriores muy separados. Las otras 2 discrepancias de esa comprobación son coincidencias de posición del mismo avión en vuelos distintos (3 h y 1,5 h de diferencia).

## 6. Repetición con 70 días (5 nov 2026)
**Resultado: 0 fallos en 61 comprobaciones.** Lo nuevo respecto al 3 nov:
- **Las secciones C y F usan solo los días de ERA5 disponibles en el entorno de Claude** (10 de 70). Prueban el código, no los datos de cada día. Los cálculos completos con los 70 días se hicieron en el Mac con el script 2.
- **Sección D, con 70 días:**
  - velocidad medida a bordo frente a distancia/tiempo: +0,85 kt (530 tramos);
  - duración del modelo frente a la observada: −0,28 % con Mach medido.
- **Sección I:** compara con el intervalo de `resultados_f3_estadistica.txt` en vez de con una cifra fija.
- **Hallazgo nuevo: el defecto del 787 depende de la dirección del viento** (`scripts/f3_direccion_b787.py`).
  - Afecta al 98–99 % de los informes del 787 con viento real del sur (150–210°), al 61–71 % con viento de 120–150° o 210–240°, y a casi ninguno en las demás direcciones.
  - Por eso su frecuencia cambia de un mes a otro (del 1 al 37 %), según cuántos días haya de viento del sur. Con los 10 días iniciales salía un 4 %; con los 70, un 13,7 %.
  - No es cosa de ciertos aviones: de 185 aviones, ninguno tiene la mitad o más de sus informes girados, y pasa igual en todos los operadores y en las tres versiones (787-8, -9 y -10).
- **Criterio del filtro ampliado:** cociente de velocidades 0,6–1,5 en vez de 0,7–1,3. Con el antiguo se escapaban 29 informes girados con viento fuerte. El RMS de ERA5 pasa de 3,44 a 3,25 m/s; la mediana no cambia (2,29). Se dice en el artículo.
- **Prueba independiente de la hora, con 70 días:** con la hora local tomada como UTC, el RMS sería de 4,96 m/s y la mediana de 3,41 m/s, frente a 3,25 y 2,29 con la hora corregida.

| Cifra | 10 días (3 nov) | **70 días (final)** |
|---|---|---|
| Tramos evaluados | 376 | **1885** |
| Exceso medio volada/óptimo | 0,30 % (0,19–0,44) | **0,23 % (0,20–0,27)** |
| Mediana | 0,13 % | **0,10 %** |
| Exceso medio por tramo | 26 kg | **22 kg** |
| Parte estructural | 0,063 % (21 %) | **0,062 % (27 %)** |
| Error vectorial de ERA5 (RMS / mediana) | 3,11 / 2,26 m/s | **3,25 / 2,29 m/s** |
| Informes del 787 con la dirección girada | 4,1 % | **13,7 %** |
