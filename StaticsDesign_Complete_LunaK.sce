// =============================================================================
//  DISEÑO ESTÁTICO COMPLETO DE EJE DE TRANSMISIÓN
//  Luna Katalina Quintero Jiménez
//  Universidad Tecnológica de Bolívar
//  Introducción al Diseño Mecánico — Prof. Edgardo Arrieta
// =============================================================================
//
//  DATOS ASIGNADOS (Hoja compartida):
//    Pot = 84 kW,  Wang = 800 RPM
//    θ1 = 16°,  θ2 = 295°
//    Diam1 = 5.0 pulg,  Diam2 = 6.5 pulg
//    φ = 23° (presión),  ψ = 19° (hélice)
//    L1 = 0.17 m,  L2 = 0.12 m,  L3 = 0.08 m
//    Sut = 305 MPa,  Num = 12
//    Presión de ajuste engranajes: 18 MPa (externa)
//
//  CRITERIOS DE FALLA:
//    a) Von Mises (Energía de distorsión)
//    b) Cortante Máximo (Teoría de Mohr)
//
// =============================================================================

clear; clc;

// ═════════════════════════════════════════════════════════════════════════════
// FUNCIONES
// ═════════════════════════════════════════════════════════════════════════════

function y = PrCruz(r, F)
    compX = r(2)*F(3) - r(3)*F(2)
    compY = r(3)*F(1) - r(1)*F(3)
    compZ = r(1)*F(2) - r(2)*F(1)
    y = [compX, compY, compZ]
endfunction

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

// Esfuerzos por presión (Lamé)
function st = Esf_Tangencial(pi_val, po_val, ri, ro, r)
    denom   = ro^2 - ri^2
    comun   = pi_val*ri^2 - po_val*ro^2
    termino = (ri*ro)^2 * (po_val - pi_val) / (r^2)
    st = (comun - termino) / denom
endfunction

function sr = Esf_Radial(pi_val, po_val, ri, ro, r)
    denom   = ro^2 - ri^2
    comun   = pi_val*ri^2 - po_val*ro^2
    termino = (ri*ro)^2 * (po_val - pi_val) / (r^2)
    sr = (comun + termino) / denom
endfunction

// Esfuerzos de flexión y torsión (eje macizo)
function sf = Esf_Flexion(Mom, diam)
    r = diam / 2
    c = r
    I = (%pi * diam^4) / 64
    sf = Mom * c / I
endfunction

function sTor = Esf_Torsion(Tor, diam)
    r = diam / 2
    I = (%pi * diam^4) / 64
    J = 2 * I
    sTor = Tor * r / J
endfunction

// Tensor completo → s1, Svm, tau_max (Mohr)
// Incluye presión externa en los engranajes (ajuste por interferencia)
function [s1, Svm, tau_max] = Calcular_Todo(Mom, Tor, diam, po_val)
    sflex = Esf_Flexion(Mom, diam)
    sTor  = Esf_Torsion(Tor, diam)

    // Presión externa: el engranaje ejerce po sobre la superficie del eje
    // Para eje macizo: ri = 0, ro = diam/2
    // En la superficie (r = ro):
    //   σ_tangencial = -po (compresión uniforme para ri=0)
    //   σ_radial     = -po
    if abs(po_val) > 0 then
        st = -po_val    // Compresión tangencial
        sr = -po_val    // Compresión radial en la superficie
    else
        st = 0
        sr = 0
    end

    // Tensor de Cauchy 3×3
    //          X(axial)  Y(radial)  Z(tangencial)
    Tensor = [sflex,       0,       sTor;
                0,        sr,         0;
              sTor,        0,        st]

    // Esfuerzos principales
    esf_ppales = real(spec(Tensor))
    esf_ppales = gsort(esf_ppales, 'g', 'd')  // s1 > s2 > s3
    S1 = esf_ppales(1)
    S2 = esf_ppales(2)
    S3 = esf_ppales(3)

    // Criterio S-normal-máx
    s1 = max(abs(esf_ppales))

    // Criterio Von Mises
    Svm = sqrt( ((S1-S2)^2 + (S2-S3)^2 + (S3-S1)^2) / 2 )

    // Criterio Cortante Máximo (Mohr)
    // τ_max = (S1 - S3) / 2
    // Falla cuando τ_max >= Sut/2
    tau_max = (S1 - S3) / 2

endfunction

// ═════════════════════════════════════════════════════════════════════════════
// DATOS DEL PROBLEMA (T00068464)
// ═════════════════════════════════════════════════════════════════════════════

mprintf("════════════════════════════════════════════════════════════════\n")
mprintf("  DISEÑO ESTÁTICO DE EJE — T00068464\n")
mprintf("  Luna Katalina Quintero Jiménez\n")
mprintf("════════════════════════════════════════════════════════════════\n\n")

// Operación
Pot_kW      = 84
n_rpm       = 800
ang_presion = 23       // φ [°]
ang_helice  = 19       // ψ [°]
pulg        = 0.0254   // Factor pulgadas → metros

// Geometría engranajes
Diam1_pulg = 5.0
Diam2_pulg = 6.5
radio1 = (Diam1_pulg * pulg) / 2    // Radio engrane 1 [m]
radio2 = (Diam2_pulg * pulg) / 2    // Radio engrane 2 [m]

// Ángulos de contacto YZ
ang_yz_1 = 16         // θ1 [°]
ang_yz_2 = 295        // θ2 [°]

// Longitudes de los tramos
L1 = 0.17    // Apoyo A → Engrane 1
L2 = 0.12    // Engrane 1 → Engrane 2
L3 = 0.08    // Engrane 2 → Apoyo B
L  = L1 + L2 + L3    // Longitud total = 0.37 m

// Posiciones axiales desde A
x1 = L1               // Engrane 1 en x = 0.17 m
x2 = L1 + L2          // Engrane 2 en x = 0.29 m

// Material
Sut_MPa = 305
Sut = Sut_MPa * 1e6   // [Pa]

// Presión de ajuste de los engranajes (externa sobre el eje)
P_ajuste = 18e6        // 18 MPa [Pa]

// Factor de seguridad deseado
Nseg_deseado = 2.5

// Imprimir datos
mprintf("--- DATOS ---\n")
mprintf("  Pot = %d kW,  n = %d RPM\n", Pot_kW, n_rpm)
mprintf("  φ = %d°,  ψ = %d°\n", ang_presion, ang_helice)
mprintf("  Diam1 = %.1f pulg (r = %.4f m)\n", Diam1_pulg, radio1)
mprintf("  Diam2 = %.1f pulg (r = %.4f m)\n", Diam2_pulg, radio2)
mprintf("  θ1 = %d°,  θ2 = %d°\n", ang_yz_1, ang_yz_2)
mprintf("  L1=%.2f, L2=%.2f, L3=%.2f → L=%.2f m\n", L1, L2, L3, L)
mprintf("  Sut = %d MPa\n", Sut_MPa)
mprintf("  Presión de ajuste = %.0f MPa\n", P_ajuste/1e6)
mprintf("  Nseg deseado = %.1f\n\n", Nseg_deseado)

// ═════════════════════════════════════════════════════════════════════════════
// ESTÁTICA: REACCIONES
// ═════════════════════════════════════════════════════════════════════════════

wang  = n_rpm * (2 * %pi / 60)
Pot_W = Pot_kW * 1000

mprintf("════════════════════════════════════════════════════════════════\n")
mprintf("  ESTÁTICA\n")
mprintf("════════════════════════════════════════════════════════════════\n")
mprintf("  ω = %.4f rad/s\n", wang)
mprintf("  Torque = %.4f N·m\n", Pot_W/wang)

// Vectores de posición
r1_ext = [x1, radio1*sind(ang_yz_1), -radio1*cosd(ang_yz_1)]
r2_ext = [x2, radio2*sind(ang_yz_2), -radio2*cosd(ang_yz_2)]

// Fuerzas: Engrane 1 entrada (+Pot), Engrane 2 salida (-Pot)
F1 = FuerzaEngrane(ang_presion, ang_helice, Pot_W, wang, radio1, ang_yz_1)
F2 = FuerzaEngrane(ang_presion, ang_helice, -Pot_W, wang, radio2, ang_yz_2)

// Momentos respecto a A
mf1 = PrCruz(r1_ext, F1)
mf2 = PrCruz(r2_ext, F2)

// Vector B
sf = -F1 - F2
sm = -mf1 - mf2

mprintf("\n  F1 = [%.2f, %.2f, %.2f] N  |F1| = %.2f N\n", F1(1), F1(2), F1(3), norm(F1))
mprintf("  F2 = [%.2f, %.2f, %.2f] N  |F2| = %.2f N\n", F2(1), F2(2), F2(3), norm(F2))
mprintf("  sf = [%.2f, %.2f, %.2f] N\n", sf(1), sf(2), sf(3))
mprintf("  sm = [%.4f, %.4f, %.4f] N·m\n", sm(1), sm(2), sm(3))

// Matriz A (6×6) — Rodillo en B (RBx = 0)
A = [1,0,0,1,0,0; 
     0,1,0,0,1,0; 
     0,0,1,0,0,1; ...
     0,0,0,0,0,0; 
     0,0,0,0,0,-L; 
     0,0,0,0,L,0]
B_vec = [sf(1); sf(2); sf(3); sm(1); sm(2); sm(3)]
A(4,4) = 1; B_vec(4) = 0;  // RBx = 0

res = A \ B_vec
RA = res(1:3)'
RB = res(4:6)'

check_F = F1 + F2 + RA + RB
mprintf("\n  RA = [%.4f, %.4f, %.4f] N  |RA| = %.4f N\n", RA(1), RA(2), RA(3), norm(RA))
mprintf("  RB = [%.4f, %.4f, %.4f] N  |RB| = %.4f N\n", RB(1), RB(2), RB(3), norm(RB))
mprintf("  ΣF = [%.6f, %.6f, %.6f]\n", check_F(1), check_F(2), check_F(3))

// ═════════════════════════════════════════════════════════════════════════════
// BARRIDO DE CARGAS INTERNAS (1000 puntos)
// ═════════════════════════════════════════════════════════════════════════════

n_cargas = 4
puntos_x = [0, x1, x2, L]
todas_F  = [RA; F1; F2; RB]
r_rel    = [[0,0,0]; r1_ext; r2_ext; [0,0,0]]
r_rel(:,1) = 0

n_puntos = 1000
dist_x   = linspace(0, L, n_puntos)
F0_mat   = zeros(n_puntos, 3)
M0_mat   = zeros(n_puntos, 3)
Mf_res   = zeros(n_puntos, 1)

for i = 1:n_puntos
    curr_x = dist_x(i)
    F_sum = [0,0,0]
    M_sum = [0,0,0]
    for j = 1:n_cargas
        if puntos_x(j) <= curr_x then
            F_sum = F_sum + todas_F(j,:)
            brazo = [puntos_x(j)-curr_x, r_rel(j,2), r_rel(j,3)]
            M_sum = M_sum + PrCruz(brazo, todas_F(j,:))
        end
    end
    F0_mat(i,:) = -F_sum
    M0_mat(i,:) = -M_sum
    Mf_res(i)   = sqrt(M0_mat(i,2)^2 + M0_mat(i,3)^2)
end

// ═════════════════════════════════════════════════════════════════════════════
// DIMENSIONAMIENTO AUTOMÁTICO (3 secciones)
// ═════════════════════════════════════════════════════════════════════════════
// Sección 1: [0, L1]       → Apoyo A a Engrane 1 (sin presión de ajuste)
// Sección 2: [L1, L1+L2]   → Engrane 1 a Engrane 2 (CON presión de ajuste)
// Sección 3: [L1+L2, L]    → Engrane 2 a Apoyo B (sin presión de ajuste)
//
// La presión de ajuste (18 MPa) actúa solo en la zona donde los
// engranajes están montados. En las secciones de los apoyos no hay
// presión externa.

mprintf("\n════════════════════════════════════════════════════════════════\n")
mprintf("  DIMENSIONAMIENTO AUTOMÁTICO\n")
mprintf("════════════════════════════════════════════════════════════════\n")

n_secciones = 3
limites_x = [0, L1, L1+L2, L]

// Presión por sección: solo la sección 2 (entre engranes) tiene presión
presion_sec = [0, P_ajuste, 0]

diametros = zeros(1, n_secciones)

for sec = 1:n_secciones
    idx_sec = find(dist_x >= limites_x(sec) & dist_x <= limites_x(sec+1))
    Mf_max  = max(Mf_res(idx_sec))
    Tor_max = max(abs(M0_mat(idx_sec, 1)))
    po_sec  = presion_sec(sec)

    mprintf("\n── Sección %d: x=[%.3f, %.3f] m ──\n", sec, limites_x(sec), limites_x(sec+1))
    mprintf("    Mf_max = %.4f N·m\n", Mf_max)
    mprintf("    Tor_max = %.4f N·m\n", Tor_max)
    mprintf("    Presión ext = %.0f MPa\n", po_sec/1e6)

    if Mf_max < 1e-10 & Tor_max < 1e-10 then
        diametros(sec) = 0.010
        mprintf("    Sin carga → D = 10 mm\n")
        continue
    end

    mprintf("\n    D[mm]  s1[MPa]  Svm[MPa]  τmax[MPa]  Nseg(VM)  Nseg(Mohr)\n")
    mprintf("    ──────────────────────────────────────────────────────────\n")

    d_prueba   = 0.005
    d_paso     = 0.001
    encontrado = %F

    while d_prueba <= 0.300
        [s1_p, Svm_p, tau_p] = Calcular_Todo(Mf_max, Tor_max, d_prueba, po_sec)

        // Nseg Von Mises: Sut / Svm
        Nvm = 99
        if Svm_p > 0 then
            Nvm = Sut / Svm_p
        end

        // Nseg Mohr (cortante máx): (Sut/2) / τ_max → Sut / (2·τ_max)
        Nmohr = 99
        if tau_p > 0 then
            Nmohr = (Sut/2) / tau_p
        end

        Ncrit = min(Nvm, Nmohr)

        if Ncrit < Nseg_deseado * 2 | Ncrit >= Nseg_deseado then
            mprintf("    %5.0f  %8.2f  %8.2f   %8.2f    %6.2f     %6.2f", ...
                    d_prueba*1000, s1_p/1e6, Svm_p/1e6, tau_p/1e6, Nvm, Nmohr)
            if Ncrit >= Nseg_deseado & ~encontrado then
                mprintf("  ← OK")
                encontrado = %T
                diametros(sec) = d_prueba
            end
            mprintf("\n")
        end

        if encontrado & d_prueba > diametros(sec) + 2*d_paso then
            break
        end
        d_prueba = d_prueba + d_paso
    end

    if ~encontrado then
        diametros(sec) = 0.300
        mprintf("    ✗ No encontrado < 300 mm\n")
    end
    mprintf("    ✓ D%d = %.0f mm\n", sec, diametros(sec)*1000)
end

// ═════════════════════════════════════════════════════════════════════════════
// RECALCULAR EN 1000 PUNTOS CON DIÁMETROS DEFINITIVOS
// ═════════════════════════════════════════════════════════════════════════════

s1_vec      = zeros(n_puntos, 1)
Svm_vec     = zeros(n_puntos, 1)
tau_max_vec = zeros(n_puntos, 1)
Nseg_vm_vec = zeros(n_puntos, 1)
Nseg_mohr_vec = zeros(n_puntos, 1)

for i = 1:n_puntos
    curr_x = dist_x(i)
    d_actual = diametros(1)
    po_actual = presion_sec(1)
    for j = 1:n_secciones
        if curr_x >= limites_x(j) & curr_x <= limites_x(j+1) then
            d_actual = diametros(j)
            po_actual = presion_sec(j)
            break
        end
    end

    Mom_local = Mf_res(i)
    Tor_local = abs(M0_mat(i, 1))

    if Mom_local > 1e-10 | Tor_local > 1e-10 then
        [s1_v, Svm_v, tau_v] = Calcular_Todo(Mom_local, Tor_local, d_actual, po_actual)
        s1_vec(i)      = s1_v
        Svm_vec(i)     = Svm_v
        tau_max_vec(i) = tau_v
    end

    if Svm_vec(i) > 0 then
        Nseg_vm_vec(i) = Sut / Svm_vec(i)
    else
        Nseg_vm_vec(i) = 99
    end
    if tau_max_vec(i) > 0 then
        Nseg_mohr_vec(i) = (Sut/2) / tau_max_vec(i)
    else
        Nseg_mohr_vec(i) = 99
    end
end

// ═════════════════════════════════════════════════════════════════════════════
// TABLA RESUMEN + VERIFICACIÓN TENSOR
// ═════════════════════════════════════════════════════════════════════════════

mprintf("\n════════════════════════════════════════════════════════════════\n")
mprintf("  TABLA DE DISEÑO FINAL\n")
mprintf("════════════════════════════════════════════════════════════════\n\n")
mprintf("  Secc  Rango [m]       D[mm]  Mf[N·m]  Tor[N·m]  P_ext  Nseg(VM)  Nseg(Mohr)\n")
mprintf("  ────────────────────────────────────────────────────────────────────────────\n")

for j = 1:n_secciones
    idx_sec = find(dist_x >= limites_x(j) & dist_x <= limites_x(j+1))
    Mf_max  = max(Mf_res(idx_sec))
    Tor_max = max(abs(M0_mat(idx_sec, 1)))

    Nvm_filt = Nseg_vm_vec(idx_sec)
    Nvm_filt(Nvm_filt >= 99) = %inf
    Nmo_filt = Nseg_mohr_vec(idx_sec)
    Nmo_filt(Nmo_filt >= 99) = %inf

    mprintf("  %3d   [%.3f,%.3f]  %5.0f  %8.2f  %8.2f  %4.0fMPa  %6.2f     %6.2f\n", ...
            j, limites_x(j), limites_x(j+1), diametros(j)*1000, ...
            Mf_max, Tor_max, presion_sec(j)/1e6, min(Nvm_filt), min(Nmo_filt))
end

// Verificación del tensor por sección
mprintf("\n════════════════════════════════════════════════════════════════\n")
mprintf("  VERIFICACIÓN: TENSOR POR SECCIÓN\n")
mprintf("════════════════════════════════════════════════════════════════\n")

for j = 1:n_secciones
    d = diametros(j)
    idx_sec = find(dist_x >= limites_x(j) & dist_x <= limites_x(j+1))
    Mf_max  = max(Mf_res(idx_sec))
    Tor_max = max(abs(M0_mat(idx_sec, 1)))
    po_sec  = presion_sec(j)

    sflex = Esf_Flexion(Mf_max, d)
    sTor  = Esf_Torsion(Tor_max, d)
    if abs(po_sec) > 0 then
        st_val = -po_sec
        sr_val = -po_sec
    else
        st_val = 0
        sr_val = 0
    end

    Tensor = [sflex, 0, sTor; 0, sr_val, 0; sTor, 0, st_val]
    esf_ppales = real(spec(Tensor))
    esf_ppales = gsort(esf_ppales, 'g', 'd')
    S1 = esf_ppales(1); S2 = esf_ppales(2); S3 = esf_ppales(3)
    Svm_val = sqrt(((S1-S2)^2+(S2-S3)^2+(S3-S1)^2)/2)
    tau_val = (S1-S3)/2

    mprintf("\n── Sección %d (D = %.0f mm) ──\n", j, d*1000)
    mprintf("    σ_flex = %.2f MPa\n", sflex/1e6)
    mprintf("    τ_tor  = %.2f MPa\n", sTor/1e6)
    if abs(po_sec) > 0 then
        mprintf("    σ_tang = %.2f MPa (presión ajuste)\n", st_val/1e6)
        mprintf("    σ_rad  = %.2f MPa (presión ajuste)\n", sr_val/1e6)
    end
    mprintf("\n    Tensor [MPa]:\n")
    mprintf("    ┌%10.2f %10.2f %10.2f┐\n", sflex/1e6, 0, sTor/1e6)
    mprintf("    │%10.2f %10.2f %10.2f│\n", 0, sr_val/1e6, 0)
    mprintf("    └%10.2f %10.2f %10.2f┘\n", sTor/1e6, 0, st_val/1e6)
    mprintf("\n    s1=%.2f, s2=%.2f, s3=%.2f MPa\n", S1/1e6, S2/1e6, S3/1e6)
    mprintf("    Svm = %.2f MPa\n", Svm_val/1e6)
    mprintf("    τ_max = %.2f MPa\n", tau_val/1e6)
    mprintf("    Nseg (Von Mises) = %.2f\n", Sut/Svm_val)
    mprintf("    Nseg (Mohr)      = %.2f\n", (Sut/2)/tau_val)
end

mprintf("\n════════════════════════════════════════════════════════════════\n")
mprintf("  DISEÑO FINAL: D1 = %.0f mm, D2 = %.0f mm, D3 = %.0f mm\n", ...
        diametros(1)*1000, diametros(2)*1000, diametros(3)*1000)
mprintf("  Material: Sut = %d MPa | Nseg = %.1f\n", Sut_MPa, Nseg_deseado)
mprintf("════════════════════════════════════════════════════════════════\n")

// ═════════════════════════════════════════════════════════════════════════════
// GRÁFICAS (10 VENTANAS)
// ═════════════════════════════════════════════════════════════════════════════

titF = ["Carga Axial (F0x)", "Cortante X-Y (F0y)", "Cortante X-Z (F0z)"]
titM = ["Torsión (M0x)", "Momento X-Z (M0y)", "Momento X-Y (M0z)"]
colF = ["k", "r", "b"]
colM = ["m", "c", "g"]

// 1-3: Fuerzas internas
for k = 1:3
    scf(k); clf();
    plot(dist_x, F0_mat(:,k), colF(k), 'linewidth', 2);
    ylabel("Fuerza [N]");
    xtitle(titF(k), "x [m]");
    xgrid(color("grey80"));
end

// 4-6: Momentos internos
for k = 1:3
    scf(k+3); clf();
    plot(dist_x, M0_mat(:,k), colM(k), 'linewidth', 2);
    ylabel("Momento [N·m]");
    xtitle(titM(k), "x [m]");
    xgrid(color("grey80"));
end

// 7: Momento flector resultante
scf(7); clf();
plot(dist_x, Mf_res, 'r-', 'linewidth', 3);
xtitle("Momento Flector Resultante", "x [m]", "Momento [N·m]");
xgrid(color("grey80"));

// 8: Esfuerzos s1 y Von Mises
scf(8); clf();
plot(dist_x, s1_vec/1e6, 'r-', 'linewidth', 2);
plot(dist_x, Svm_vec/1e6, 'b-', 'linewidth', 2);
Sperm = Sut / Nseg_deseado;
plot([0, L], [Sperm/1e6, Sperm/1e6], 'k--', 'linewidth', 1);
xstring(L*0.5, Sperm/1e6*1.05, sprintf("Sut/%.1f = %.0f MPa", Nseg_deseado, Sperm/1e6))
legend(["s1", "Svm (Von Mises)", "Esf. permisible"], 'in_upper_left');
xtitle("Esfuerzos", "x [m]", "[MPa]");
xgrid(color("grey80"));

// 9: Factor de seguridad (Von Mises y Mohr)
scf(9); clf();
Nvm_plot  = min(Nseg_vm_vec, 10);
Nmo_plot  = min(Nseg_mohr_vec, 10);
plot(dist_x, Nvm_plot, 'b-', 'linewidth', 2);
plot(dist_x, Nmo_plot, 'r-', 'linewidth', 2);
plot([0, L], [Nseg_deseado, Nseg_deseado], 'k--', 'linewidth', 1.5);
xstring(L*0.01, Nseg_deseado+0.15, sprintf("Nseg = %.1f", Nseg_deseado))
plot([0, L], [2, 2], 'g--', 'linewidth', 1);
xstring(L*0.7, 2.15, "Nseg = 2")
legend(["Von Mises", "Mohr (τ_max)", "Nseg deseado"], 'in_upper_left');
xtitle("Factor de Seguridad", "x [m]", "Nseg");
xgrid(color("grey80"));

// 10: Perfil del eje
scf(10); clf();
for j = 1:n_secciones
    r_mm = diametros(j)/2 * 1000;
    xi = limites_x(j);
    xf = limites_x(j+1);
    plot([xi, xf], [r_mm, r_mm], 'b-', 'linewidth', 3);
    plot([xi, xf], [-r_mm, -r_mm], 'b-', 'linewidth', 3);
    if j > 1 then
        r_prev = diametros(j-1)/2 * 1000;
        plot([xi, xi], [r_prev, r_mm], 'b-', 'linewidth', 2);
        plot([xi, xi], [-r_prev, -r_mm], 'b-', 'linewidth', 2);
    end
    xstring((xi+xf)/2, r_mm+1, sprintf("D%d=%dmm", j, diametros(j)*1000))
end
plot([0, 0], [-diametros(1)/2*1000, diametros(1)/2*1000], 'b-', 'linewidth', 2);
plot([L, L], [-diametros(3)/2*1000, diametros(3)/2*1000], 'b-', 'linewidth', 2);
plot([0, L], [0, 0], 'k--', 'linewidth', 0.5);
// Marcar posiciones de engranajes
plot([x1, x1], [-diametros(2)/2*1000-3, diametros(2)/2*1000+3], 'g--', 'linewidth', 1);
plot([x2, x2], [-diametros(2)/2*1000-3, diametros(2)/2*1000+3], 'g--', 'linewidth', 1);
xstring(x1, diametros(2)/2*1000+4, "Eng.1")
xstring(x2, diametros(2)/2*1000+4, "Eng.2")
xtitle("Perfil del Eje Diseñado", "x [m]", "Radio [mm]");
xgrid(color("grey80"));

mprintf("\n>> 10 ventanas generadas.\n")
mprintf(">> 1-3: Fuerzas | 4-6: Momentos | 7: Mf_res\n")
mprintf(">> 8: Esfuerzos | 9: Nseg (VM + Mohr) | 10: Perfil\n")
