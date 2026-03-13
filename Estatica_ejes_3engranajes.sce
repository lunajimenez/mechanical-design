// =============================================================================
// DISEÑO ESTÁTICO DE EJES: 3 ENGRANAJES
// Autora: Luna Katalina Quintero Jiménez
// Universidad Tecnológica de Bolívar
// Metodología del Prof. Edgardo Arrieta
// =============================================================================
// Configuración de potencia:
//   Engrane 1: +0.4·P (entrada parcial)
//   Engrane 2: -P     (salida total)
//   Engrane 3: +0.6·P (entrada parcial)
//
// Balance de torque: 0.4P + 0.6P - P = 0  ✓
// =============================================================================

// En este caso, lo que cambia es la cantidad de fuentes de carga y cómo se reparte la potencia.

clear; clc;

// =====================================================================
// SECCIÓN 1: ENTRADA DE DATOS
// =====================================================================

mprintf("================================================================\n")
mprintf("  DISEÑO ESTÁTICO DE EJES: 3 ENGRANAJES\n")
mprintf("  Distribución: 0.4P (entrada) / -P (salida) / 0.6P (entrada)\n")
mprintf("================================================================\n\n")

Pot_kW      = input("  Potencia total P [kW]      : ")
n           = input("  Velocidad de rotación [RPM]: ")
ang_presion = input("  Ángulo de presión φ [°]    : ")
ang_helice  = input("  Ángulo de hélice ψ [°]     : ")
L           = input("  Longitud total A→B [m]     : ")

mprintf("\n--- ENGRANE 1 (entrada 0.4P) ---\n")
radio1   = input("  Radio primitivo [m]        : ")
ang_yz_1 = input("  Ángulo de contacto YZ [°]  : ")
x1       = input("  Posición axial desde A [m] : ")

mprintf("\n--- ENGRANE 2 (salida -P) ---\n")
radio2   = input("  Radio primitivo [m]        : ")
ang_yz_2 = input("  Ángulo de contacto YZ [°]  : ")
x2       = input("  Posición axial desde A [m] : ")

mprintf("\n--- ENGRANE 3 (entrada 0.6P) ---\n")
radio3   = input("  Radio primitivo [m]        : ")
ang_yz_3 = input("  Ángulo de contacto YZ [°]  : ")
x3       = input("  Posición axial desde A [m] : ")

rodillo_en = input("\n  ¿Rodillo en B o A? (A/B): ", "string")

// =====================================================================
// SECCIÓN 2: FUNCIONES
// =====================================================================

function fu = FuerzaEngrane(ang_presion, ang_helice, Pot, wang, radio, ang_yz)
    Torque = Pot / wang
    Ft = Torque / radio
    Fa = Ft * tand(ang_helice)
    Fr = Ft * tand(ang_presion)

    Ur = [0, cosd(ang_yz), sind(ang_yz)]
    Ut = [0, -sind(ang_yz), cosd(ang_yz)]
    Ua = [1, 0, 0]

    fu = -Fr*Ur + Ft*Ut + Fa*Ua
endfunction

function y = PrCruz(r, F)
    compX = r(2)*F(3) - r(3)*F(2)
    compY = r(3)*F(1) - r(1)*F(3)
    compZ = r(1)*F(2) - r(2)*F(1)
    y = [compX, compY, compZ]
endfunction

// =====================================================================
// SECCIÓN 3: CÁLCULO DE REACCIONES
// =====================================================================

wang  = n * (2 * %pi / 60)
Pot_W = Pot_kW * 1000

// --- Vectores de posición desde A (origen) ---
r1_ext = [x1, radio1*sind(ang_yz_1), -radio1*cosd(ang_yz_1)]
r2_ext = [x2, radio2*sind(ang_yz_2), -radio2*cosd(ang_yz_2)]
r3_ext = [x3, radio3*sind(ang_yz_3), -radio3*cosd(ang_yz_3)]

// --- Fuerzas en cada engranaje ---
// Engrane 1: entrada parcial (+0.4P)
// Engrane 2: salida total   (-P)
// Engrane 3: entrada parcial (+0.6P)
F1 = FuerzaEngrane(ang_presion, ang_helice,  0.4*Pot_W, wang, radio1, ang_yz_1)
F2 = FuerzaEngrane(ang_presion, ang_helice, -1.0*Pot_W, wang, radio2, ang_yz_2)
F3 = FuerzaEngrane(ang_presion, ang_helice,  0.6*Pot_W, wang, radio3, ang_yz_3)

// --- Momentos respecto a A ---
mf1 = PrCruz(r1_ext, F1)
mf2 = PrCruz(r2_ext, F2)
mf3 = PrCruz(r3_ext, F3)

// --- Vector B (lado derecho) ---
// ΣF = 0 → RA + RB = -F1 - F2 - F3
// ΣM = 0 → análogamente
sf = -F1 - F2 - F3
sm = -mf1 - mf2 - mf3

// --- Matriz A (6×6) ---
//             RAx  RAy  RAz  RBx  RBy   RBz
A = [           1,   0,   0,   1,   0,    0;
                0,   1,   0,   0,   1,    0;
                0,   0,   1,   0,   0,    1;
                0,   0,   0,   0,   0,    0;
                0,   0,   0,   0,   0,   -L;
                0,   0,   0,   0,   L,    0 ]

B = [sf(1); sf(2); sf(3); sm(1); sm(2); sm(3)]

// --- Condición de rodillo ---
if convstr(rodillo_en, "u") == "B" then
    A(4,4) = 1; B(4) = 0;
else
    A(4,1) = 1; B(4) = 0;
end

// --- Resolución ---
res = A \ B
RA = res(1:3)'
RB = res(4:6)'

// --- Verificación ---
check_F = F1(:)' + F2(:)' + F3(:)' + RA(:)' + RB(:)'
mprintf("\n>> Verificación ΣF (debe ≈ 0): [%.6f, %.6f, %.6f]\n", ...
        check_F(1), check_F(2), check_F(3))

// --- Presentación de reacciones ---
mprintf("\n--- FUERZAS EN LOS ENGRANAJES ---\n")
mprintf("  F1 (0.4P) = [%.4f, %.4f, %.4f] N  |F1| = %.4f N\n", F1(1), F1(2), F1(3), norm(F1))
mprintf("  F2 (-P)   = [%.4f, %.4f, %.4f] N  |F2| = %.4f N\n", F2(1), F2(2), F2(3), norm(F2))
mprintf("  F3 (0.6P) = [%.4f, %.4f, %.4f] N  |F3| = %.4f N\n", F3(1), F3(2), F3(3), norm(F3))

// Explicación del uso de: norm. 
// norm calcula la magnitud (o módulo) del vector
// Lo usamos para mostrar la fuerza total que siente el eje en ese punto, como un solo número en newtons 
// Las tres componentes dicen hacia dónde empuja la fuerza, pero la magnitud dice cuánto empuja en total 
// norm no es parte del cálculo ni afecta ningún resultado

mprintf("\n--- REACCIONES ---\n")
mprintf("  RA = [%.4f, %.4f, %.4f] N  |RA| = %.4f N\n", RA(1), RA(2), RA(3), norm(RA))
mprintf("  RB = [%.4f, %.4f, %.4f] N  |RB| = %.4f N\n", RB(1), RB(2), RB(3), norm(RB))

// --- Balance de torque ---
mprintf("\n--- BALANCE DE TORQUE ---\n")
mprintf("  sm(1) = %.6f N·m (debe ≈ 0)\n", sm(1))

// =====================================================================
// SECCIÓN 4: CÁLCULO DE DIAGRAMAS (BARRIDO CON 5 CARGAS)
// =====================================================================
// Ahora hay 5 puntos de carga: A, Engrane 1, Engrane 2, Engrane 3, B

puntos_x = [0, x1, x2, x3, L]

todas_F = [RA; F1; F2; F3; RB]

r_rel = [[0,0,0]; r1_ext; r2_ext; r3_ext; [0,0,0]]
r_rel(:,1) = 0

n_puntos = 1000
dist_x   = linspace(0, L, n_puntos)

F0_mat = zeros(n_puntos, 3)
M0_mat = zeros(n_puntos, 3)
Mf_res = zeros(n_puntos, 1)

for i = 1:n_puntos
    curr_x = dist_x(i)
    F_sum_left = [0, 0, 0]
    M_sum_left = [0, 0, 0]

    // Ahora se revisan 5 cargas en vez de 4
    for j = 1:5
        if puntos_x(j) <= curr_x then
            F_sum_left = F_sum_left + todas_F(j, :)
            brazo = [puntos_x(j) - curr_x, r_rel(j, 2), r_rel(j, 3)]
            M_sum_left = M_sum_left + PrCruz(brazo, todas_F(j, :))
        end
    end

    F0_mat(i, :) = -F_sum_left
    M0_mat(i, :) = -M_sum_left
    Mf_res(i) = sqrt(M0_mat(i, 2)^2 + M0_mat(i, 3)^2)
end

// =====================================================================
// SECCIÓN 5: GENERACIÓN DE GRÁFICAS
// =====================================================================

titulos = ["Carga Axial (F0x)", "Cortante X-Y (F0y)", "Cortante X-Z (F0z)", ..
           "Torsión (M0x)", "Momento X-Z (M0y)", "Momento X-Y (M0z)", ..
           "Momento Flector Resultante"]

colores = ["k", "r", "b", "m", "c", "g", "r"]

for k = 1:6
    scf(k);
    clf();

    if k <= 3 then
        plot(dist_x, F0_mat(:,k), colores(k), 'linewidth', 2);
        ylabel("Fuerza [N]");
    else
        plot(dist_x, M0_mat(:,k-3), colores(k), 'linewidth', 2);
        ylabel("Momento [N·m]");
    end

    xtitle(titulos(k), "Posición axial x [m]");
    xgrid(color("grey80"));
end

scf(7);
clf();
plot(dist_x, Mf_res, 'r-', 'linewidth', 3);
xtitle(titulos(7), "Posición axial x [m]", "Momento [N·m]");
xgrid(color("grey80"));

mprintf("\n>> Proceso completado. Revisa las 7 ventanas de gráficos.\n")
