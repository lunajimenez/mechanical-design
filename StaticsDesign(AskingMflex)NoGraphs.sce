// =============================================================================
//  DIMENSIONAMIENTO DE EJES POR SECCIONES
//  Entrada directa: Mf, Tor, Sut, Nseg → Salida: Diámetro mínimo
//  Autora: Luna Katalina Quintero Jiménez y equipo
//  Universidad Tecnológica de Bolívar
//  Metodología del Prof. Edgardo Arrieta
// =============================================================================
//
//  USO:
//    El usuario ingresa para cada sección del eje:
//      - Momento flector resultante [N·m]
//      - Torque (torsión) [N·m]
//      - Sut del material [MPa]
//      - Nseg deseado
//
//    El código calcula el diámetro mínimo que cumple ambos criterios:
//      - S-normal-máx:  Nseg = Sut / s1
//      - Von Mises:     Nseg = Sut / Svm
//
//    Soporta ejes macizos (espesor=0) o huecos (espesor>0),
//    con o sin presión interna/externa.
//
// =============================================================================

// Datos importantes para una primera prueba: 
// pfalla = 1e-4
// Nseg = 2.23
// Sut = 340MPa
// Mflex1 = 640, Torsión = 0
// Mflex2 = 580, Torsión = 750
// Mflex1 = 1200, Torsión = 0

// =============================================================================

clear; clc;

// ═════════════════════════════════════════════════════════════════════════════
// SECCIÓN 1: FUNCIONES
// ═════════════════════════════════════════════════════════════════════════════

// --- Esfuerzo tangencial (Lamé, cilindro a presión) ---
function st = Esf_Tangencial(pi_val, po_val, ri, ro, r)
    denom   = ro^2 - ri^2
    comun   = pi_val*ri^2 - po_val*ro^2
    termino = (ri*ro)^2 * (po_val - pi_val) / (r^2)
    st = (comun - termino) / denom
endfunction

// --- Esfuerzo radial (Lamé, cilindro a presión) ---
function sr = Esf_Radial(pi_val, po_val, ri, ro, r)
    denom   = ro^2 - ri^2
    comun   = pi_val*ri^2 - po_val*ro^2
    termino = (ri*ro)^2 * (po_val - pi_val) / (r^2)
    sr = (comun + termino) / denom
endfunction

// --- Esfuerzo de flexión ---
function sf = Esf_Flexion(Mom, diam, espesor)
    r = diam / 2
    if espesor > 0 then
        ri = r - espesor
    else
        ri = 0
    end
    c = r
    I = (%pi*r^4)/4 - (%pi*ri^4)/4
    sf = Mom * c / I
endfunction

// --- Esfuerzo de torsión ---
function sTor = Esf_Torsion(Tor, diam, espesor)
    r = diam / 2
    if espesor > 0 then
        ri = r - espesor
    else
        ri = 0
    end
    I = (%pi*r^4)/4 - (%pi*ri^4)/4
    J = 2 * I
    sTor = Tor * r / J
endfunction

// --- Tensor completo → s1 y Svm ---
//
// TENSOR DE CAUCHY:
//     ┌ σ_flex    0    τ_tor ┐
//     │   0      σ_r     0   │
//     └ τ_tor     0    σ_t   ┘
//
// Sin presión (pi=po=0): σ_r = σ_t = 0 → tensor simplificado de clase
//
function [s1, Svm] = Calcular_Esfuerzos(Mom, Tor, diam, espesor, pi_val, po_val)

    sflex = Esf_Flexion(Mom, diam, espesor)
    sTor  = Esf_Torsion(Tor, diam, espesor)

    r  = diam / 2
    if espesor > 0 then
        ri = r - espesor
    else
        ri = 0
    end
    ro = r

    if (abs(pi_val) > 0 | abs(po_val) > 0) & espesor > 0 then
        st = Esf_Tangencial(pi_val, po_val, ri, ro, ro)
        sr = Esf_Radial(pi_val, po_val, ri, ro, ro)
    else
        st = 0
        sr = 0
    end

    Tensor = [sflex,  0,    sTor;
                0,    sr,     0;
              sTor,    0,    st]

    esf_ppales = real(spec(Tensor))
    esf_ppales = gsort(esf_ppales, 'g', 'd')
    S1 = esf_ppales(1)
    S2 = esf_ppales(2)
    S3 = esf_ppales(3)

    s1  = max(abs(esf_ppales))
    Svm = sqrt( ((S1-S2)^2 + (S2-S3)^2 + (S3-S1)^2) / 2 )

endfunction

// ═════════════════════════════════════════════════════════════════════════════
// SECCIÓN 2: ENTRADA DE DATOS
// ═════════════════════════════════════════════════════════════════════════════

mprintf("════════════════════════════════════════════════════════════════\n")
mprintf("  DIMENSIONAMIENTO DE EJES POR SECCIONES\n")
mprintf("  Entrada: Mf, Tor, Sut, Nseg → Salida: Diámetro mínimo\n")
mprintf("════════════════════════════════════════════════════════════════\n\n")

// --- Tipo de eje ---
mprintf("--- TIPO DE EJE ---\n")
mprintf("  Macizo: espesor = 0\n")
mprintf("  Hueco:  espesor > 0 (ej: 0.00635 m = 1/4 pulg)\n")
espesor = input("  Espesor de pared [m] (0 si macizo): ")

// --- Presión (solo si aplica) ---
mprintf("\n--- PRESIÓN (0 si no aplica) ---\n")
pi_val = input("  Presión interna [Pa]: ")
po_val = input("  Presión externa [Pa]: ")

// --- Número de secciones ---
mprintf("\n--- SECCIONES DEL EJE ---\n")
n_secciones = input("  Número de secciones a dimensionar: ")

// --- Datos por sección ---
Mom_sec     = zeros(1, n_secciones)
Tor_sec     = zeros(1, n_secciones)
Sut_sec     = zeros(1, n_secciones)
Nseg_sec    = zeros(1, n_secciones)

for i = 1:n_secciones
    mprintf("\n══ SECCIÓN %d ══\n", i)
    Mom_sec(i)  = input("  Momento flector resultante Mf [N·m]: ")
    Tor_sec(i)  = input("  Torque (torsión) [N·m]             : ")
    Sut_sec(i)  = input("  Sut del material [MPa]             : ") * 1e6  // → Pa
    Nseg_sec(i) = input("  Factor de seguridad deseado         : ")
end

// ═════════════════════════════════════════════════════════════════════════════
// SECCIÓN 3: CÁLCULO DE DIÁMETROS
// ═════════════════════════════════════════════════════════════════════════════
//
// Para cada sección se itera el diámetro desde 5 mm en pasos de 1 mm
// hasta encontrar el mínimo que cumple:
//   Nseg(s1)  = Sut / s1  >= Nseg_deseado
//   Nseg(Svm) = Sut / Svm >= Nseg_deseado
//
// Esto replica la iteración manual del profesor (Diap. 28-30):
//   D=15 → Nseg=1.26 (no), D=18 → Nseg=2.28 (no), D=19 → Nseg=2.56 (sí)
//
// ═════════════════════════════════════════════════════════════════════════════

mprintf("\n\n════════════════════════════════════════════════════════════════\n")
mprintf("  RESULTADOS\n")
mprintf("════════════════════════════════════════════════════════════════\n")

diametros = zeros(1, n_secciones)

for i = 1:n_secciones

    mprintf("\n── Sección %d ──\n", i)
    mprintf("    Mf   = %.2f N·m\n", Mom_sec(i))
    mprintf("    Tor  = %.2f N·m\n", Tor_sec(i))
    mprintf("    Sut  = %.0f MPa\n", Sut_sec(i)/1e6)
    mprintf("    Nseg deseado = %.1f\n", Nseg_sec(i))

    // Caso especial: sin carga → diámetro mínimo constructivo
    if Mom_sec(i) < 1e-10 & Tor_sec(i) < 1e-10 then
        mprintf("    Sin carga → D mínimo constructivo = 10 mm\n")
        diametros(i) = 0.010
        continue
    end

    // --- Iteración ---
    d_prueba    = 0.005     // Empieza en 5 mm
    d_paso      = 0.001     // Sube de 1 mm en 1 mm
    d_max       = 0.500     // Límite 500 mm
    encontrado  = %F

    // Mostrar tabla de iteraciones (como hace el profesor)
    mprintf("\n    D [mm]   s1 [MPa]   Svm [MPa]   Nseg(s1)   Nseg(VM)\n")
    mprintf("    ─────────────────────────────────────────────────────\n")

    while d_prueba <= d_max

        // Verificar espesor válido
        if espesor > 0 & espesor >= d_prueba/2 then
            d_prueba = d_prueba + d_paso
            continue
        end

        [s1_val, Svm_val] = Calcular_Esfuerzos( ...
            Mom_sec(i), Tor_sec(i), d_prueba, espesor, pi_val, po_val)

        if s1_val > 0 then
            Nseg_s1 = Sut_sec(i) / s1_val
        else
            Nseg_s1 = 99
        end

        if Svm_val > 0 then
            Nseg_vm = Sut_sec(i) / Svm_val
        else
            Nseg_vm = 99
        end

        Nseg_critico = min(Nseg_s1, Nseg_vm)

        // Imprimir solo los diámetros cercanos a la solución
        // (para no llenar la consola con 500 líneas)
        if Nseg_critico < Nseg_sec(i) * 2 | Nseg_critico >= Nseg_sec(i) then
            mprintf("    %5.0f    %8.2f    %8.2f     %6.2f     %6.2f", ...
                    d_prueba*1000, s1_val/1e6, Svm_val/1e6, Nseg_s1, Nseg_vm)
            if Nseg_critico >= Nseg_sec(i) & ~encontrado then
                mprintf("  ← CUMPLE")
                encontrado = %T
                diametros(i) = d_prueba
            end
            mprintf("\n")
        end

        // Seguir mostrando un par más después de encontrar
        if encontrado then
            // Mostrar 2 diámetros más para comparación
            if d_prueba > diametros(i) + 2*d_paso then
                break
            end
        end

        d_prueba = d_prueba + d_paso
    end

    if ~encontrado then
        diametros(i) = d_max
        mprintf("    ✗ No se encontró D < %.0f mm\n", d_max*1000)
    end

    mprintf("\n    ✓ DIÁMETRO SELECCIONADO: D = %.0f mm\n", diametros(i)*1000)
end

// ═════════════════════════════════════════════════════════════════════════════
// SECCIÓN 4: TABLA RESUMEN
// ═════════════════════════════════════════════════════════════════════════════

mprintf("\n\n════════════════════════════════════════════════════════════════\n")
mprintf("  TABLA DE DISEÑO FINAL\n")
mprintf("════════════════════════════════════════════════════════════════\n\n")
mprintf("  Sección   Mf [N·m]   Tor [N·m]   Sut [MPa]   Nseg    D [mm]\n")
mprintf("  ────────────────────────────────────────────────────────────\n")

for i = 1:n_secciones
    mprintf("  %4d      %8.2f    %8.2f     %6.0f     %4.1f     %5.0f\n", ...
            i, Mom_sec(i), Tor_sec(i), Sut_sec(i)/1e6, Nseg_sec(i), diametros(i)*1000)
end

mprintf("  ────────────────────────────────────────────────────────────\n")
if espesor > 0 then
    mprintf("  Eje hueco, espesor = %.2f mm\n", espesor*1000)
else
    mprintf("  Eje macizo\n")
end
if pi_val > 0 | po_val > 0 then
    mprintf("  Presión interna: %.2f MPa, externa: %.2f MPa\n", pi_val/1e6, po_val/1e6)
end

mprintf("\n════════════════════════════════════════════════════════════════\n")

// ═════════════════════════════════════════════════════════════════════════════
// SECCIÓN 5: VERIFICACIÓN DETALLADA (Tensor y esfuerzos principales)
// ═════════════════════════════════════════════════════════════════════════════
// Imprime el tensor de Cauchy y los esfuerzos principales para cada
// sección con el diámetro seleccionado, para incluir en el informe.
// ═════════════════════════════════════════════════════════════════════════════

mprintf("\n\n════════════════════════════════════════════════════════════════\n")
mprintf("  VERIFICACIÓN DETALLADA POR SECCIÓN\n")
mprintf("════════════════════════════════════════════════════════════════\n")

for i = 1:n_secciones

    d = diametros(i)
    r = d / 2

    if espesor > 0 then
        ri = r - espesor
    else
        ri = 0
    end

    // Recalcular esfuerzos con diámetro final
    sflex = Esf_Flexion(Mom_sec(i), d, espesor)
    sTor  = Esf_Torsion(Tor_sec(i), d, espesor)

    if (abs(pi_val) > 0 | abs(po_val) > 0) & espesor > 0 then
        st = Esf_Tangencial(pi_val, po_val, ri, r, r)
        sr = Esf_Radial(pi_val, po_val, ri, r, r)
    else
        st = 0
        sr = 0
    end

    Tensor = [sflex, 0, sTor; 0, sr, 0; sTor, 0, st]

    esf_ppales = real(spec(Tensor))
    esf_ppales = gsort(esf_ppales, 'g', 'd')
    S1 = esf_ppales(1)
    S2 = esf_ppales(2)
    S3 = esf_ppales(3)
    s1_final  = max(abs(esf_ppales))
    Svm_final = sqrt(((S1-S2)^2 + (S2-S3)^2 + (S3-S1)^2)/2)

    Nseg_s1_final = Sut_sec(i) / s1_final
    Nseg_vm_final = Sut_sec(i) / Svm_final

    mprintf("\n── Sección %d (D = %.0f mm) ──\n", i, d*1000)
    mprintf("    I = %.4e m⁴\n", (%pi*r^4)/4 - (%pi*ri^4)/4)
    mprintf("    J = %.4e m⁴\n", 2*((%pi*r^4)/4 - (%pi*ri^4)/4))
    mprintf("    σ_flex = %.2f MPa\n", sflex/1e6)
    mprintf("    τ_tor  = %.2f MPa\n", sTor/1e6)
    if st ~= 0 | sr ~= 0 then
        mprintf("    σ_tang = %.2f MPa\n", st/1e6)
        mprintf("    σ_rad  = %.2f MPa\n", sr/1e6)
    end

    mprintf("\n    Tensor de Cauchy [MPa]:\n")
    mprintf("    ┌ %10.2f  %10.2f  %10.2f ┐\n", sflex/1e6, 0, sTor/1e6)
    mprintf("    │ %10.2f  %10.2f  %10.2f │\n", 0, sr/1e6, 0)
    mprintf("    └ %10.2f  %10.2f  %10.2f ┘\n", sTor/1e6, 0, st/1e6)

    mprintf("\n    Esfuerzos principales [MPa]:\n")
    mprintf("    s1 = %.2f,  s2 = %.2f,  s3 = %.2f\n", S1/1e6, S2/1e6, S3/1e6)
    mprintf("    |s1| máx = %.2f MPa\n", s1_final/1e6)
    mprintf("    Svm      = %.2f MPa\n", Svm_final/1e6)

    mprintf("\n    Nseg (S-n-máx)   = %.2f\n", Nseg_s1_final)
    mprintf("    Nseg (Von Mises) = %.2f\n", Nseg_vm_final)
end

mprintf("\n════════════════════════════════════════════════════════════════\n")
mprintf("  Proceso completado.\n")
mprintf("════════════════════════════════════════════════════════════════\n")
