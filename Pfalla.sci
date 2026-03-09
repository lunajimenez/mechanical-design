// --- CÁLCULO DE PROBABILIDAD DE FALLA ---

// 1. Datos de la Resistencia del Material (SUT)
mu_SUT = 320;            
CX_SUT = 0.08;           
desv_SUT = CX_SUT * mu_SUT; 

// 2. Datos del Esfuerzo (calculados en Excel/NX)
mu_esfuerzo = 65.7215;       
CX_esfuerzo = 0.1286;       
desv_esfuerzo = CX_esfuerzo * mu_esfuerzo; 

// 3. Variable "Resta" (Margen de Seguridad: Resistencia - Esfuerzo)
mu_R = mu_SUT - mu_esfuerzo; 
desv_R = sqrt(desv_esfuerzo^2 + desv_SUT^2); 

// 4. Cálculo de la Probabilidad de Falla
// Evaluamos en 0 porque si (Resistencia - Esfuerzo) <= 0, el componente falla.
P = cdfnor("PQ", 0, mu_R, desv_R);

// 5. Resultados
printf("La probabilidad de falla es: %f (o %.2f%%)\n", P, P*100);
