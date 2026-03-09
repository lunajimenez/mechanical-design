// =============================================================================
// DISEÑO ESTÁTICO DE EJES: CÁLCULO DE REACCIONES Y DIAGRAMAS
// Autora: Luna Katalina Quintero Jiménez
// Universidad Tecnológica de Bolívar
// Metodología del Prof. Edgardo Arrieta

// =============================================================================
// Proceso global (Diapositiva 3):
//   DCL → Reacciones → Diag. Cortante → Diag. Momentos → Mom. Resultante
//
// Este script cubre las primeras cuatro etapas de ese proceso.
// =============================================================================

clear; clc;

// =====================================================================
// SECCIÓN 1: ENTRADA DE DATOS
// =====================================================================
// Se solicitan al usuario todos los parámetros necesarios para que
// el mismo script resuelva cualquier configuración de eje sin
// modificar el código fuente.
// =====================================================================

mprintf("================================================================\n")
mprintf("  DISEÑO ESTÁTICO DE EJES: CÁLCULO Y DIAGRAMAS\n")
mprintf("================================================================\n\n")

Pot_kW      = input("  Potencia [kW]              : ")
n           = input("  Velocidad de rotación [RPM]: ")
ang_presion = input("  Ángulo de presión φ [°]    : ")
ang_helice  = input("  Ángulo de hélice ψ [°]     : ")
L           = input("  Longitud total A→B [m]     : ")

mprintf("\n--- ENGRANE 1 ---\n")
radio1   = input("  Radio primitivo [m]        : ")
ang_yz_1 = input("  Ángulo de contacto YZ [°]  : ")
x1       = input("  Posición axial desde A [m] : ")

mprintf("\n--- ENGRANE 2 ---\n")
radio2   = input("  Radio primitivo [m]        : ")
ang_yz_2 = input("  Ángulo de contacto YZ [°]  : ")
x2       = input("  Posición axial desde A [m] : ")

rodillo_en = input("\n  ¿Rodillo en B o A? (A/B): ", "string")

// =====================================================================
// SECCIÓN 2: FUNCIONES
// =====================================================================
// Ambas funciones retornan vectores FILA [1x3] de forma explícita.
// Esto es crítico porque el operador A\B de Scilab retorna vectores
// COLUMNA [3x1], y si se mezclan dimensiones [1x3] + [3x1] el código
// produce un error o resultados incorrectos de [3x3].
// =====================================================================

// --- FuerzaEngrane ---------------------------------------------------
// Convierte potencia + geometría del engranaje en un vector de fuerza 3D.
//
// Referencia: Diapositiva 6
//   Los vectores unitarios Ur, Ut, Ua dependen del ángulo de contacto
//   θ en el plano YZ. Esto permite rotar la fuerza automáticamente
//   a la posición correcta sin proyecciones manuales.
//
// Signo de Pot:
//   Pot > 0 → el engranaje ENTREGA potencia al eje (entrada)
//   Pot < 0 → el eje ENTREGA potencia al engranaje (salida)
//   El cambio de signo invierte la dirección de Ft y Fa.
// -----------------------------------------------------------------
function fu = FuerzaEngrane(ang_presion, ang_helice, Pot, wang, radio, ang_yz)
    Torque = Pot / wang
    Ft = Torque / radio            // Fuerza tangencial [N]
    Fa = Ft * tand(ang_helice)     // Fuerza axial [N]
    Fr = Ft * tand(ang_presion)    // Fuerza radial [N]

    // Vectores unitarios en el punto de contacto (Diap. 6)
    Ur = [0, cosd(ang_yz), sind(ang_yz)]    // Radial: del eje hacia el diente
    Ut = [0, -sind(ang_yz), cosd(ang_yz)]   // Tangencial: perpendicular a Ur en YZ
    Ua = [1, 0, 0]                          // Axial: a lo largo del eje (+X)

    // Resultado: vector fila [1x3]
    // -Fr*Ur: la fuerza radial empuja HACIA el centro del eje
    // +Ft*Ut: la tangencial actúa en la dirección de giro
    // +Fa*Ua: la axial empuja a lo largo del eje
    fu = -Fr*Ur + Ft*Ut + Fa*Ua
endfunction

// --- PrCruz ----------------------------------------------------------
// Producto cruz manual: y = r × F
// Retorna EXPLÍCITAMENTE un vector fila [1x3] para evitar el error
// de incompatibilidad dimensional [1x3] + [3x1].
//
// Componentes (regla de la mano derecha):
//   X: torque (rotación alrededor del eje)
//   Y: flexión en el plano XZ
//   Z: flexión en el plano XY
// -----------------------------------------------------------------
function y = PrCruz(r, F)
    compX = r(2)*F(3) - r(3)*F(2)
    compY = r(3)*F(1) - r(1)*F(3)
    compZ = r(1)*F(2) - r(2)*F(1)
    y = [compX, compY, compZ]    // Fila explícita [1x3]
endfunction

// =====================================================================
// SECCIÓN 3: CÁLCULO DE REACCIONES
// =====================================================================
// Referencia: Diapositivas 7-14
//
// Se construye el sistema Ax = B donde:
//   x = [RAx, RAy, RAz, RBx, RBy, RBz]'
//   Filas 1-3: ΣFx, ΣFy, ΣFz = 0     (Diap. 8)
//   Filas 5-6: ΣMy, ΣMz = 0           (Diap. 9)
//   Fila 4:    condición de rodillo    (Diap. 14)
//
// Momentos respecto a A (origen):
//   rA = [0,0,0] → no genera momento
//   rB = [L,0,0] → rB × RB = [0, -L·RBz, L·RBy]  (Diap. 11)
// =====================================================================

wang  = n * (2 * %pi / 60)    // RPM → rad/s
Pot_W = Pot_kW * 1000         // kW → W

// --- Vectores de posición desde A (origen) hasta punto de contacto ---
// Convención: r = [distancia_axial, radio·sin(θ), -radio·cos(θ)]
// La componente Z lleva signo negativo por la convención del DCL
r1_ext = [x1, radio1*sind(ang_yz_1), -radio1*cosd(ang_yz_1)]
r2_ext = [x2, radio2*sind(ang_yz_2), -radio2*cosd(ang_yz_2)]

// --- Fuerzas en cada engranaje ---
F1 = FuerzaEngrane(ang_presion, ang_helice,  Pot_W, wang, radio1, ang_yz_1)  // Entrada
F2 = FuerzaEngrane(ang_presion, ang_helice, -Pot_W, wang, radio2, ang_yz_2)  // Salida

// --- Momentos de las fuerzas externas respecto a A ---
mf1 = PrCruz(r1_ext, F1)
mf2 = PrCruz(r2_ext, F2)

// --- Vector B (lado derecho) ---
// De ΣF = 0: RA + RB + F1 + F2 = 0 → RA + RB = -F1 - F2
// De ΣM = 0: análogamente para momentos
sf = -F1 - F2        // Sumatoria de fuerzas conocidas [N]
sm = -mf1 - mf2      // Sumatoria de momentos conocidos [N·m]

// --- Matriz A (6×6) ---
// Referencia: Diapositivas 12-14
//             RAx  RAy  RAz  RBx  RBy   RBz
A = [           1,   0,   0,   1,   0,    0;   // ΣFx = 0
                0,   1,   0,   0,   1,    0;   // ΣFy = 0
                0,   0,   1,   0,   0,    1;   // ΣFz = 0
                0,   0,   0,   0,   0,    0;   // ΣMx → degenerada (Diap. 13: fila de ceros)
                0,   0,   0,   0,   0,   -L;   // ΣMy = 0 → coeficiente -L en RBz
                0,   0,   0,   0,   L,    0 ]  // ΣMz = 0 → coeficiente +L en RBy

B = [sf(1); sf(2); sf(3); sm(1); sm(2); sm(3)]

// --- Condición de rodillo (reemplaza la fila 4 degenerada) ---
// Diapositiva 14: se fija RBx=0 (rodillo en B) o RAx=0 (rodillo en A)
// para que el sistema no sea singular.
if convstr(rodillo_en, "u") == "B" then
    A(4,4) = 1; B(4) = 0;    // RBx = 0
else
    A(4,1) = 1; B(4) = 0;    // RAx = 0
end

// --- Resolución ---
res = A \ B         // Eliminación gaussiana (Scilab nativo)
RA = res(1:3)'      // Transpuesta a fila [1x3] — crítico para compatibilidad
RB = res(4:6)'      // Transpuesta a fila [1x3]

// --- Verificación rápida ---
check_F = F1(:)' + F2(:)' + RA(:)' + RB(:)'
mprintf("\n>> Verificación ΣF (debe ≈ 0): [%.6f, %.6f, %.6f]\n", ...
        check_F(1), check_F(2), check_F(3))

// --- Presentación de reacciones ---
mprintf("\n  RA = [%.4f, %.4f, %.4f] N\n", RA(1), RA(2), RA(3))
mprintf("  RB = [%.4f, %.4f, %.4f] N\n", RB(1), RB(2), RB(3))

// =====================================================================
// SECCIÓN 4: CÁLCULO DE DIAGRAMAS (BARRIDO DE CARGAS INTERNAS)
// =====================================================================
//
// FUNDAMENTO TEÓRICO (Diapositiva 16):
// ─────────────────────────────────────
// Imagina que cortas el eje con una sierra en la posición x.
// Para que la mitad izquierda no se caiga ni gire, la cara del
// corte debe ejercer una fuerza interna F0 y un momento interno M0.
//
//   F0(x) = −Σ (todas las fuerzas a la izquierda de x)
//   M0(x) = −Σ (todos los momentos de esas fuerzas respecto a x)
//
// El signo negativo aparece porque F0 y M0 son las fuerzas y momentos
// que la mitad DERECHA ejerce sobre la mitad IZQUIERDA para mantener
// el equilibrio. Si la suma de fuerzas a la izquierda empuja hacia
// arriba, F0 debe empujar hacia abajo → signo opuesto.
//
// METODOLOGÍA NUMÉRICA (Diapositiva 17):
// ──────────────────────────────────────
// En vez de resolver analíticamente por tramos (que requiere escribir
// ecuaciones diferentes entre cada par de cargas), se usa un BARRIDO
// NUMÉRICO: se evalúan F0 y M0 en 1000 puntos uniformemente
// distribuidos a lo largo del eje. En cada punto se pregunta:
// "¿cuáles cargas están a mi izquierda?" y se suman sus contribuciones.
//
// CARGAS QUE ACTÚAN SOBRE EL EJE:
// ────────────────────────────────
// Hay exactamente 4 puntos de carga: Apoyo A, Engrane 1, Engrane 2, Apoyo B.
// Cada uno tiene:
//   - Una posición axial (puntos_x)
//   - Un vector de fuerza 3D (todas_F)
//   - Un vector de brazo radial (r_rel): la distancia en Y,Z desde el
//     centro del eje hasta donde actúa la fuerza
//
// Para los apoyos, r_rel = [0,0,0] porque las reacciones actúan
// directamente sobre la línea central del eje.
// Para los engranes, r_rel tiene componentes Y,Z que representan
// el radio primitivo proyectado según el ángulo de contacto.
//
// =====================================================================

// --- Posiciones axiales de las 4 cargas ---
puntos_x = [0, x1, x2, L]

// --- Fuerzas en cada punto de carga ---
// Fila 1: Reacción en A (actúa en x=0)
// Fila 2: Fuerza del engrane 1 (actúa en x=x1)
// Fila 3: Fuerza del engrane 2 (actúa en x=x2)
// Fila 4: Reacción en B (actúa en x=L)
todas_F = [RA; F1; F2; RB]

// --- Brazos radiales (distancia en Y,Z desde el centro del eje) ---
// Para los apoyos: [0,0,0] (la fuerza actúa en el centro)
// Para los engranes: las componentes Y,Z del vector de posición
//   (la componente X se maneja aparte con puntos_x)
r_rel = [[0,0,0]; r1_ext; r2_ext; [0,0,0]]
r_rel(:,1) = 0   // Anulamos la componente X porque la posición axial
                  // ya está en puntos_x. Si no se anulara, el brazo
                  // axial se contaría dos veces.

// --- Malla de evaluación ---
n_puntos = 1000                      // Resolución del barrido
dist_x   = linspace(0, L, n_puntos)  // 1000 puntos de 0 a L

// --- Matrices de almacenamiento ---
// F0_mat: fuerza interna en cada punto [Vx, Vy, Vz] (N)
//   Col 1 (Vx): Carga axial interna
//   Col 2 (Vy): Fuerza cortante en el plano XY
//   Col 3 (Vz): Fuerza cortante en el plano XZ
F0_mat = zeros(n_puntos, 3)

// M0_mat: momento interno en cada punto [Mx, My, Mz] (N·m)
//   Col 1 (Mx): Torsión (torque que el eje transmite en esa sección)
//   Col 2 (My): Momento flector en el plano XZ
//   Col 3 (Mz): Momento flector en el plano XY
M0_mat = zeros(n_puntos, 3)

// Mf_res: momento flector resultante (escalar)
//   Es la magnitud total de la flexión: sqrt(My² + Mz²)
//   Este valor se usa para calcular el esfuerzo flector σ = Mf·c/I
Mf_res = zeros(n_puntos, 1)

// ─── INICIO DEL BARRIDO ─────────────────────────────────────────────
// Para cada punto x_i del eje (1000 iteraciones):
for i = 1:n_puntos
    curr_x = dist_x(i)            // Posición del corte imaginario
    F_sum_left = [0, 0, 0]        // Acumulador de fuerzas a la izquierda
    M_sum_left = [0, 0, 0]        // Acumulador de momentos a la izquierda

    // Revisamos las 4 cargas (A, Engrane 1, Engrane 2, B)
    for j = 1:4
        if puntos_x(j) <= curr_x then
            // ──────────────────────────────────────────────────────
            // Esta carga está a la IZQUIERDA del corte → contribuye
            // ──────────────────────────────────────────────────────

            // 1) Sumar la fuerza
            F_sum_left = F_sum_left + todas_F(j, :)

            // 2) Calcular el BRAZO DE PALANCA desde la carga hasta el corte
            //
            //    brazo_x = posición_carga - posición_corte (NEGATIVO,
            //              porque la carga está a la izquierda)
            //    brazo_y = r_rel(j,2) = componente Y del radio de contacto
            //    brazo_z = r_rel(j,3) = componente Z del radio de contacto
            //
            //    NOTA: el brazo se mide DESDE la carga HACIA el corte,
            //    por eso brazo_x es negativo. Esto es consistente con
            //    la convención de la Diapositiva 16: "el brazo de palanca
            //    se mide desde el punto de aplicación de la carga externa
            //    hasta la cara del corte imaginario."
            brazo = [puntos_x(j) - curr_x, r_rel(j, 2), r_rel(j, 3)]

            // 3) Momento = brazo × fuerza (producto cruz)
            //    Cada componente del resultado tiene significado físico:
            //    - Comp X: contribución al TORQUE
            //      (fuerzas tangenciales × radio del engranaje)
            //    - Comp Y: contribución a la FLEXIÓN en el plano XZ
            //      (fuerzas en Z × distancia axial al corte)
            //    - Comp Z: contribución a la FLEXIÓN en el plano XY
            //      (fuerzas en Y × distancia axial al corte)
            M_sum_left = M_sum_left + PrCruz(brazo, todas_F(j, :))
        end
    end

    // ──────────────────────────────────────────────────────────────
    // APLICAR EL SIGNO NEGATIVO (Diapositiva 16):
    // "F0(x) es el NEGATIVO de la suma de todas las fuerzas
    //  a la izquierda del punto x"
    //
    // ¿POR QUÉ negativo?
    // Porque F0 es la fuerza que la mitad DERECHA del eje ejerce
    // sobre la mitad IZQUIERDA. Si a la izquierda hay una fuerza
    // neta hacia arriba, la cara del corte debe empujar hacia ABAJO
    // para mantener el equilibrio → signo contrario.
    // ──────────────────────────────────────────────────────────────
    F0_mat(i, :) = -F_sum_left
    M0_mat(i, :) = -M_sum_left

    // ──────────────────────────────────────────────────────────────
    // MOMENTO FLECTOR RESULTANTE (Diapositiva 17):
    // Mf = sqrt(My² + Mz²)
    //
    // Combina la flexión de ambos planos en un solo escalar.
    // Este es el valor que se usa en la ecuación de esfuerzo flector:
    //   σ_f = Mf · c / I
    // donde c es la distancia al centroide e I el momento de inercia.
    // La sección más crítica del eje será donde Mf_res sea máximo.
    // ──────────────────────────────────────────────────────────────
    Mf_res(i) = sqrt(M0_mat(i, 2)^2 + M0_mat(i, 3)^2)
end
// ─── FIN DEL BARRIDO ─────────────────────────────────────────────────

// =====================================================================
// SECCIÓN 5: GENERACIÓN DE GRÁFICAS
// =====================================================================
// Referencia: Diapositiva 17
// Se generan 7 ventanas independientes (una por cada diagrama).
// Se usan ventanas separadas porque las escalas de Newton [N] y
// Newton-metro [N·m] son muy diferentes y solaparlas confundiría
// la lectura.
//
// Gráficas 1-3: Fuerzas internas (F0)
//   1. Carga Axial F0x     → fuerza normal a la sección transversal
//   2. Cortante XY  F0y    → fuerza cortante en el plano XY
//   3. Cortante XZ  F0z    → fuerza cortante en el plano XZ
//
// Gráficas 4-6: Momentos internos (M0)
//   4. Torsión M0x         → torque que transmite el eje
//   5. Momento XZ  M0y     → flexión en el plano XZ
//   6. Momento XY  M0z     → flexión en el plano XY
//
// Gráfica 7: Momento Flector Resultante
//   sqrt(M0y² + M0z²)      → magnitud total de flexión
//   Esta es la gráfica más importante para el diseño del eje.
// =====================================================================

titulos = ["Carga Axial (F0x)", "Cortante X-Y (F0y)", "Cortante X-Z (F0z)", ..
           "Torsión (M0x)", "Momento X-Z (M0y)", "Momento X-Y (M0z)", ..
           "Momento Flector Resultante"]

colores = ["k", "r", "b", "m", "c", "g", "r"]

for k = 1:6
    scf(k);
    clf();

    if k <= 3 then
        // Gráficas de FUERZA interna (columnas 1, 2, 3 de F0_mat)
        plot(dist_x, F0_mat(:,k), colores(k), 'linewidth', 2);
        ylabel("Fuerza [N]");
    else
        // Gráficas de MOMENTO interno (columnas 1, 2, 3 de M0_mat)
        plot(dist_x, M0_mat(:,k-3), colores(k), 'linewidth', 2);
        ylabel("Momento [N·m]");
    end

    xtitle(titulos(k), "Posición axial x [m]");
    xgrid(color("grey80"));
end

// --- Ventana 7: Momento Flector Resultante ---
scf(7);
clf();
plot(dist_x, Mf_res, 'r-', 'linewidth', 3);
xtitle(titulos(7), "Posición axial x [m]", "Momento [N·m]");
xgrid(color("grey80"));

// --- Resultado final ---
mprintf("\n>> Proceso completado. Revisa las 7 ventanas de gráficos.\n")

// =====================================================================
// RESUMEN DE LO QUE CALCULA CADA GRÁFICA:
// =====================================================================
//
// GRÁFICA 1 - Carga Axial (F0x):
//   Fuerza interna a lo largo del eje. Cambia solo donde hay
//   componentes axiales (Fa de los engranes o reacciones axiales).
//   Produce esfuerzo normal: σ_a = F0x / A
//
// GRÁFICA 2 - Cortante XY (F0y):
//   Fuerza cortante en el plano XY. Salta en cada carga que tenga
//   componente Y. Produce esfuerzo cortante: τ = V·Q / (I·b)
//
// GRÁFICA 3 - Cortante XZ (F0z):
//   Igual que la anterior pero en el plano XZ.
//
// GRÁFICA 4 - Torsión (M0x):
//   Torque interno que transmite el eje. Es constante entre los dos
//   engranes (donde se transmite potencia) y cero fuera de ellos.
//   Produce esfuerzo cortante: τ = T·c / J
//
// GRÁFICA 5 - Momento XZ (M0y):
//   Momento flector que dobla el eje en el plano XZ.
//   Varía linealmente entre cargas puntuales.
//
// GRÁFICA 6 - Momento XY (M0z):
//   Momento flector que dobla el eje en el plano XY.
//
// GRÁFICA 7 - Momento Flector Resultante:
//   sqrt(M0y² + M0z²). Combina ambas flexiones.
//   El MÁXIMO de esta curva indica la sección crítica del eje.
//   Se usa para: σ_f = Mf_max · c / I → criterio de Von Mises.
//
// =====================================================================
