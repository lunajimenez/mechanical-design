// =============================================================================
//  PARCIAL 2A — DISEÑO ESTÁTICO DE EJE CON PROBABILIDAD DE FALLA
//  T00068464 — Luna Katalina Quintero Jiménez
//  Universidad Tecnológica de Bolívar — Dismec 202610
//  Prof. Edgardo Arrieta
// =============================================================================
//
//  DATOS ASIGNADOS:
//    Pot=55kW, Wang=850RPM, θ1=29°, θ2=190°
//    Diam1=5.9", Diam2=6.5", φ=23°, ψ=19°
//    L1=0.17m, L2=0.12m, L3=0.08m, Sut=233MPa
//    Cx_sut=0.065, Cx_esfuerzos=0.12, dados por el profesor
//    PFALLA = 10^(-6), debido a mi código T00068464
//    Criterio: Von Mises
//    Para este código necesitamos el Tensor COMPLETO: incluye cortantes Txy, Txz
//
//  TEORÍA DE PROBABILIDAD DE FALLA:
//  ─────────────────────────────────
//  La resistencia S y el esfuerzo σ son variables aleatorias normales.
//
//    μ_S = Sut (media de la resistencia)
//    σ_S = Cx_sut × μ_S (desviación estándar de la resistencia)
//
//    μ_σ = Svm (media del esfuerzo Von Mises calculado)
//    σ_σ = Cx_esf × μ_σ (desviación estándar del esfuerzo)
//
//  La variable Z = S - σ también es normal con:
//    μ_Z = μ_S - μ_σ
//    σ_Z = sqrt(σ_S² + σ_σ²)
//
//  La falla ocurre cuando Z < 0. La probabilidad de falla es:
//  PFALLA = 10^(-6) -> z = 4.7534
//  Lo determinamos haciendo: z = -cdfnor("X", 0, 1, 1e-6, 1-1e-6) en consola.
//
//  CONDICIÓN DE DISEÑO:
//    μ_Z / σ_Z >= z_requerido
//    (μ_S - μ_σ) / sqrt(σ_S² + σ_σ²) >= 4.7534
//
//  Esto reemplaza al Nseg.
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

function [Svm, s1] = Tensor_Completo(Mf, Tor, Vy, Vz, diam, po_val)
    r = diam / 2
    A = %pi * r^2
    I = (%pi * diam^4) / 64
    J = 2 * I
    c = r

    if Mf > 1e-10 then
        sflex = Mf * c / I
    else
        sflex = 0
    end

    if abs(Tor) > 1e-10 then
        tau_tor = Tor * r / J
    else
        tau_tor = 0
    end

    tau_xy = (4/3) * abs(Vy) / A
    tau_xz_V = (4/3) * abs(Vz) / A
    tau_xz = tau_tor + tau_xz_V

    if abs(po_val) > 0 then
        st = -po_val
        sr = -po_val
    else
        st = 0
        sr = 0
    end

    Tensor = [sflex,      tau_xy,      tau_xz;
              tau_xy,       sr,           0;
              tau_xz,        0,          st]

    esf_ppales = real(spec(Tensor))
    esf_ppales = gsort(esf_ppales, 'g', 'd')
    S1 = esf_ppales(1)
    S2 = esf_ppales(2)
    S3 = esf_ppales(3)

    s1 = max(abs(esf_ppales))
    Svm = sqrt( ((S1-S2)^2 + (S2-S3)^2 + (S3-S1)^2) / 2 )
endfunction

function [z_val, Pfalla_ok] = Verificar_Pfalla(Svm, mu_S, Cx_sut, Cx_esf, z_req)
    mu_sigma = Svm
    sigma_S  = Cx_sut * mu_S
    sigma_sigma = Cx_esf * mu_sigma

    mu_Z    = mu_S - mu_sigma
    sigma_Z = sqrt(sigma_S^2 + sigma_sigma^2)

    if sigma_Z > 0 then
        z_val = mu_Z / sigma_Z
    else
        z_val = 99
    end

    Pfalla_ok = (z_val >= z_req)
endfunction

// ═════════════════════════════════════════════════════════════════════════════
// DATOS DEL PROBLEMA (T00068464)
// ═════════════════════════════════════════════════════════════════════════════

mprintf("================================================================\n")
mprintf("  PARCIAL 2A - DISENO CON PROBABILIDAD DE FALLA\n")
mprintf("  T00068464 - Luna Katalina Quintero Jimenez\n")
mprintf("================================================================\n\n")

Pot_kW      = 55
n_rpm       = 850
ang_presion = 23
ang_helice  = 19
pulg        = 0.0254

radio1 = (5.9 * pulg) / 2
radio2 = (6.5 * pulg) / 2
ang_yz_1 = 29
ang_yz_2 = 190

L1 = 0.17;  L2 = 0.12;  L3 = 0.08
L  = L1 + L2 + L3
x1 = L1
x2 = L1 + L2

Sut_MPa = 233
Sut = Sut_MPa * 1e6
Cx_sut = 0.065
Cx_esf = 0.12
Pfalla_objetivo = 1e-6

z_requerido = 4.7534

P_ajuste = 18e6

mprintf("--- DATOS ---\n")
mprintf("  Pot=%d kW, n=%d RPM, fi=%d, psi=%d\n", Pot_kW, n_rpm, ang_presion, ang_helice)
mprintf("  R1=%.4f m (5.9 pulg), R2=%.4f m (6.5 pulg)\n", radio1, radio2)
mprintf("  theta1=%d, theta2=%d\n", ang_yz_1, ang_yz_2)
mprintf("  L1=%.2f, L2=%.2f, L3=%.2f -> L=%.2f m\n", L1, L2, L3, L)
mprintf("  Sut=%d MPa, Cx_sut=%.3f, Cx_esf=%.2f\n", Sut_MPa, Cx_sut, Cx_esf)
mprintf("  PFALLA objetivo = 10^(-6) -> z_req = %.4f\n", z_requerido)
mprintf("  Presion ajuste = %.0f MPa\n\n", P_ajuste/1e6)

// ═════════════════════════════════════════════════════════════════════════════
// ESTATICA
// ═════════════════════════════════════════════════════════════════════════════

wang  = n_rpm * (2 * %pi / 60)
Pot_W = Pot_kW * 1000

mprintf("================================================================\n")
mprintf("  ESTATICA\n")
mprintf("================================================================\n")
mprintf("  w = %.4f rad/s,  T = %.4f N.m\n", wang, Pot_W/wang)

r1_ext = [x1, radio1*sind(ang_yz_1), -radio1*cosd(ang_yz_1)]
r2_ext = [x2, radio2*sind(ang_yz_2), -radio2*cosd(ang_yz_2)]

F1 = FuerzaEngrane(ang_presion, ang_helice, Pot_W, wang, radio1, ang_yz_1)
F2 = FuerzaEngrane(ang_presion, ang_helice, -Pot_W, wang, radio2, ang_yz_2)

mf1 = PrCruz(r1_ext, F1)
mf2 = PrCruz(r2_ext, F2)

sf = -F1 - F2
sm = -mf1 - mf2

A_mat = [1,0,0,1,0,0; 0,1,0,0,1,0; 0,0,1,0,0,1; ...
         0,0,0,0,0,0; 0,0,0,0,0,-L; 0,0,0,0,L,0]
B_vec = [sf(1); sf(2); sf(3); sm(1); sm(2); sm(3)]
A_mat(4,4) = 1; B_vec(4) = 0;

res = A_mat \ B_vec
RA = res(1:3)'
RB = res(4:6)'

check_F = F1 + F2 + RA + RB
mprintf("\n  F1 = [%.2f, %.2f, %.2f] N\n", F1(1), F1(2), F1(3))
mprintf("  F2 = [%.2f, %.2f, %.2f] N\n", F2(1), F2(2), F2(3))
mprintf("  RA = [%.4f, %.4f, %.4f] N\n", RA(1), RA(2), RA(3))
mprintf("  RB = [%.4f, %.4f, %.4f] N\n", RB(1), RB(2), RB(3))
mprintf("  sumF = [%.6f, %.6f, %.6f]\n", check_F(1), check_F(2), check_F(3))

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
// DIMENSIONAMIENTO CON PROBABILIDAD DE FALLA
// ═════════════════════════════════════════════════════════════════════════════

mprintf("\n================================================================\n")
mprintf("  DIMENSIONAMIENTO (PFALLA <= 10^(-6), z_req = %.4f)\n", z_requerido)
mprintf("================================================================\n")

n_secciones = 3
limites_x = [0, L1, L1+L2, L]
presion_sec = [P_ajuste, P_ajuste, 0]

diametros = zeros(1, n_secciones)

for sec = 1:n_secciones
    idx_sec = find(dist_x >= limites_x(sec) & dist_x <= limites_x(sec+1))
    Mf_max  = max(Mf_res(idx_sec))
    Tor_max = max(abs(M0_mat(idx_sec, 1)))
    Vy_max  = max(abs(F0_mat(idx_sec, 2)))
    Vz_max  = max(abs(F0_mat(idx_sec, 3)))
    po_sec  = presion_sec(sec)

    mprintf("\n-- Seccion %d: x=[%.3f, %.3f] m --\n", sec, limites_x(sec), limites_x(sec+1))
    mprintf("    Mf_max  = %.4f N.m\n", Mf_max)
    mprintf("    Tor_max = %.4f N.m\n", Tor_max)
    mprintf("    Vy_max  = %.4f N\n", Vy_max)
    mprintf("    Vz_max  = %.4f N\n", Vz_max)
    mprintf("    P_ext   = %.0f MPa\n", po_sec/1e6)

    if Mf_max < 1e-10 & Tor_max < 1e-10 & Vy_max < 1e-10 & Vz_max < 1e-10 then
        diametros(sec) = 0.010
        mprintf("    Sin carga -> D = 10 mm\n")
        continue
    end

    mprintf("\n    D[mm]  Svm[MPa]  z_calc   z_req    Estado\n")
    mprintf("    ----------------------------------------------\n")

    d_prueba   = 0.005
    d_paso     = 0.001
    encontrado = %F

    while d_prueba <= 0.300
        [Svm_p, s1_p] = Tensor_Completo(Mf_max, Tor_max, Vy_max, Vz_max, d_prueba, po_sec)

        [z_calc, ok] = Verificar_Pfalla(Svm_p, Sut, Cx_sut, Cx_esf, z_requerido)

        if Svm_p > 0 then
            Nseg_eq = Sut / Svm_p
        else
            Nseg_eq = 99
        end

        if z_calc < z_requerido * 2 | z_calc >= z_requerido then
            mprintf("    %5.0f  %8.2f  %6.2f   %6.2f", ...
                    d_prueba*1000, Svm_p/1e6, z_calc, z_requerido)
            if ok & ~encontrado then
                mprintf("   <- CUMPLE (Nseg_eq=%.2f)", Nseg_eq)
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
        mprintf("    X No encontrado < 300 mm\n")
    end
    mprintf("    >> D%d = %.0f mm\n", sec, diametros(sec)*1000)
end

// ═════════════════════════════════════════════════════════════════════════════
// RECALCULAR EN 1000 PUNTOS CON DIAMETROS DEFINITIVOS
// ═════════════════════════════════════════════════════════════════════════════

Svm_vec     = zeros(n_puntos, 1)
s1_vec      = zeros(n_puntos, 1)
z_vec       = zeros(n_puntos, 1)
Nseg_eq_vec = zeros(n_puntos, 1)

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

    Mf_i  = Mf_res(i)
    Tor_i = abs(M0_mat(i, 1))
    Vy_i  = abs(F0_mat(i, 2))
    Vz_i  = abs(F0_mat(i, 3))

    if Mf_i > 1e-10 | Tor_i > 1e-10 | Vy_i > 1e-10 | Vz_i > 1e-10 then
        [Svm_v, s1_v] = Tensor_Completo(Mf_i, Tor_i, Vy_i, Vz_i, d_actual, po_actual)
        Svm_vec(i) = Svm_v
        s1_vec(i)  = s1_v
        [z_v, dummy] = Verificar_Pfalla(Svm_v, Sut, Cx_sut, Cx_esf, z_requerido)
        z_vec(i) = z_v
    else
        z_vec(i) = 99
    end

    if Svm_vec(i) > 0 then
        Nseg_eq_vec(i) = Sut / Svm_vec(i)
    else
        Nseg_eq_vec(i) = 99
    end
end

// ═════════════════════════════════════════════════════════════════════════════
// TABLA RESUMEN Y VERIFICACION DEL TENSOR
// ═════════════════════════════════════════════════════════════════════════════

mprintf("\n================================================================\n")
mprintf("  TABLA DE DISENO FINAL\n")
mprintf("================================================================\n\n")
mprintf("  Secc  D[mm]  Mf[N.m]  Tor[N.m]  Vy[N]    Vz[N]    Svm[MPa]  z_calc\n")
mprintf("  --------------------------------------------------------------------\n")

for j = 1:n_secciones
    idx_sec = find(dist_x >= limites_x(j) & dist_x <= limites_x(j+1))
    Mf_max  = max(Mf_res(idx_sec))
    Tor_max = max(abs(M0_mat(idx_sec, 1)))
    Vy_max  = max(abs(F0_mat(idx_sec, 2)))
    Vz_max  = max(abs(F0_mat(idx_sec, 3)))

    [Svm_final, dummy] = Tensor_Completo(Mf_max, Tor_max, Vy_max, Vz_max, diametros(j), presion_sec(j))
    [z_final, dummy] = Verificar_Pfalla(Svm_final, Sut, Cx_sut, Cx_esf, z_requerido)

    mprintf("  %3d   %5.0f  %8.2f  %8.2f  %7.2f  %7.2f  %8.2f   %5.2f\n", ...
            j, diametros(j)*1000, Mf_max, Tor_max, Vy_max, Vz_max, Svm_final/1e6, z_final)
end

// Verificacion detallada del tensor
mprintf("\n================================================================\n")
mprintf("  VERIFICACION: TENSOR COMPLETO POR SECCION\n")
mprintf("================================================================\n")

for j = 1:n_secciones
    d = diametros(j)
    r = d/2
    A_sec = %pi * r^2
    I_sec = (%pi * d^4) / 64
    J_sec = 2 * I_sec

    idx_sec = find(dist_x >= limites_x(j) & dist_x <= limites_x(j+1))
    Mf_max  = max(Mf_res(idx_sec))
    Tor_max = max(abs(M0_mat(idx_sec, 1)))
    Vy_max  = max(abs(F0_mat(idx_sec, 2)))
    Vz_max  = max(abs(F0_mat(idx_sec, 3)))
    po_sec  = presion_sec(j)

    sflex    = Mf_max * r / I_sec
    tau_tor  = Tor_max * r / J_sec
    tau_xy   = (4/3) * Vy_max / A_sec
    tau_xz_V = (4/3) * Vz_max / A_sec
    tau_xz   = tau_tor + tau_xz_V

    if abs(po_sec) > 0 then
        st = -po_sec; sr = -po_sec
    else
        st = 0; sr = 0
    end

    Tensor = [sflex, tau_xy, tau_xz; tau_xy, sr, 0; tau_xz, 0, st]
    esf_ppales = real(spec(Tensor))
    esf_ppales = gsort(esf_ppales, 'g', 'd')
    S1 = esf_ppales(1); S2 = esf_ppales(2); S3 = esf_ppales(3)
    Svm_v = sqrt(((S1-S2)^2+(S2-S3)^2+(S3-S1)^2)/2)

    sigma_S = Cx_sut * Sut
    sigma_sigma = Cx_esf * Svm_v
    mu_Z = Sut - Svm_v
    sigma_Z = sqrt(sigma_S^2 + sigma_sigma^2)
    z_val = mu_Z / sigma_Z

    mprintf("\n-- Seccion %d (D = %.0f mm) --\n", j, d*1000)
    mprintf("    I = %.4e m4,  J = %.4e m4,  A = %.4e m2\n", I_sec, J_sec, A_sec)
    mprintf("    s_flex   = %.2f MPa\n", sflex/1e6)
    mprintf("    t_tor    = %.2f MPa\n", tau_tor/1e6)
    mprintf("    t_xy(V)  = %.2f MPa  (cortante por Vy)\n", tau_xy/1e6)
    mprintf("    t_xz(V)  = %.2f MPa  (cortante por Vz)\n", tau_xz_V/1e6)
    mprintf("    t_xz_tot = %.2f MPa  (torsion + cortante Vz)\n", tau_xz/1e6)
    if abs(po_sec) > 0 then
        mprintf("    s_r = s_t = %.2f MPa  (presion ajuste)\n", sr/1e6)
    end

    mprintf("\n    Tensor [MPa]:\n")
    mprintf("    |%10.2f %10.2f %10.2f|\n", sflex/1e6, tau_xy/1e6, tau_xz/1e6)
    mprintf("    |%10.2f %10.2f %10.2f|\n", tau_xy/1e6, sr/1e6, 0)
    mprintf("    |%10.2f %10.2f %10.2f|\n", tau_xz/1e6, 0, st/1e6)

    mprintf("\n    s1=%.2f, s2=%.2f, s3=%.2f MPa\n", S1/1e6, S2/1e6, S3/1e6)
    mprintf("    Svm = %.2f MPa\n", Svm_v/1e6)
    mprintf("    Nseg_eq = Sut/Svm = %.2f\n", Sut/Svm_v)

    mprintf("\n    --- Probabilidad de falla ---\n")
    mprintf("    mu_S = %.0f MPa,  sigma_S = %.2f MPa\n", Sut/1e6, sigma_S/1e6)
    mprintf("    mu_sigma = %.2f MPa,  sigma_sigma = %.2f MPa\n", Svm_v/1e6, sigma_sigma/1e6)
    mprintf("    mu_Z = %.2f MPa,  sigma_Z = %.2f MPa\n", mu_Z/1e6, sigma_Z/1e6)
    mprintf("    z = mu_Z/sigma_Z = %.4f\n", z_val)
    mprintf("    z_req = %.4f  ->  ", z_requerido)
    if z_val >= z_requerido then
        mprintf("PFALLA <= 10^(-6) OK\n")
    else
        mprintf("NO CUMPLE\n")
    end
end

mprintf("\n================================================================\n")
mprintf("  DISENO FINAL: D1=%.0fmm, D2=%.0fmm, D3=%.0fmm\n", ...
        diametros(1)*1000, diametros(2)*1000, diametros(3)*1000)
mprintf("  Criterio: Von Mises + Pfalla <= 10^(-6)\n")
mprintf("  Tensor: COMPLETO (incluye t_xy, t_xz por cortantes)\n")
mprintf("================================================================\n")

// ═════════════════════════════════════════════════════════════════════════════
// GRAFICAS (10 VENTANAS)
// ═════════════════════════════════════════════════════════════════════════════

titF = ["Carga Axial (F0x)", "Cortante X-Y (F0y)", "Cortante X-Z (F0z)"]
titM = ["Torsion (M0x)", "Momento X-Z (M0y)", "Momento X-Y (M0z)"]
colF = ["k", "r", "b"]
colM = ["m", "c", "g"]

for k = 1:3
    scf(k); clf();
    plot(dist_x, F0_mat(:,k), colF(k), 'linewidth', 2);
    ylabel("Fuerza [N]"); xtitle(titF(k), "x [m]"); xgrid(color("grey80"));
end

for k = 1:3
    scf(k+3); clf();
    plot(dist_x, M0_mat(:,k), colM(k), 'linewidth', 2);
    ylabel("Momento [N.m]"); xtitle(titM(k), "x [m]"); xgrid(color("grey80"));
end

scf(7); clf();
plot(dist_x, Mf_res, 'r-', 'linewidth', 3);
xtitle("Momento Flector Resultante", "x [m]", "Momento [N.m]");
xgrid(color("grey80"));

scf(8); clf();
plot(dist_x, Svm_vec/1e6, 'b-', 'linewidth', 2);
plot(dist_x, s1_vec/1e6, 'r-', 'linewidth', 1.5);
xtitle("Esfuerzos: Von Mises y s1", "x [m]", "[MPa]");
legend(["Svm (Von Mises)", "s1 (ppal max)"], 'in_upper_left');
xgrid(color("grey80"));

scf(9); clf();
z_plot = min(z_vec, 15);
plot(dist_x, z_plot, 'b-', 'linewidth', 2);
plot([0, L], [z_requerido, z_requerido], 'r--', 'linewidth', 1.5);
xstring(L*0.02, z_requerido+0.3, "z_req = 4.75 (Pfalla=1e-6)")
xtitle("Indice de confiabilidad z", "x [m]", "z = mu_Z / sigma_Z");
xgrid(color("grey80"));

scf(10); clf();
for j = 1:n_secciones
    r_mm = diametros(j)/2 * 1000;
    xi = limites_x(j); xf = limites_x(j+1);
    plot([xi, xf], [r_mm, r_mm], 'b-', 'linewidth', 3);
    plot([xi, xf], [-r_mm, -r_mm], 'b-', 'linewidth', 3);
    if j > 1 then
        r_prev = diametros(j-1)/2 * 1000;
        plot([xi, xi], [r_prev, r_mm], 'b-', 'linewidth', 2);
        plot([xi, xi], [-r_prev, -r_mm], 'b-', 'linewidth', 2);
    end
    xstring((xi+xf)/2, r_mm+1, "D"+string(j)+"="+string(diametros(j)*1000)+"mm")
end
plot([0, 0], [-diametros(1)/2*1000, diametros(1)/2*1000], 'b-', 'linewidth', 2);
plot([L, L], [-diametros(3)/2*1000, diametros(3)/2*1000], 'b-', 'linewidth', 2);
plot([0, L], [0, 0], 'k--', 'linewidth', 0.5);
xtitle("Perfil del Eje", "x [m]", "Radio [mm]");
xgrid(color("grey80"));

mprintf("\n>> 10 ventanas generadas.\n")
