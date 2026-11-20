# =====================================================================
# 3o TRABALHO PRATICO - HIDROLOGIA ESTOCASTICA (EHR806)
# QUESTAO 11 - AR(2) com phi11 = 0,9 e rho2 = 0,8. Determinar phi22.
# =====================================================================
graphics.off()

# --- TESTE PRIMEIRO com o exemplo do proprio slide: r1=0,57 e r2=0,07
r1 <- 0.57; r2 <- 0.07
cat("Exemplo do slide: phi12 =", round(r1*(1-r2)/(1-r1^2), 3),
    " phi22 =", round((r2 - r1^2)/(1 - r1^2), 3), " (slide: 0,785 e -0,38)\n")

# --- A QUESTAO ----------------------------------------------------------
phi11 <- 0.9
rho2  <- 0.8
rho1  <- phi11            # para o AR(1), phi11 = rho1  (esta nos slides)

# Formula do slide
phi22 <- (rho2 - rho1^2)/(1 - rho1^2)
phi12 <- rho1*(1 - rho2)/(1 - rho1^2)
cat("phi22 =", round(phi22, 4), "\n")
cat("phi12 =", round(phi12, 4), "\n")

# Conferencia 1: pela forma matricial do slide (Phi = inversa(P) %*% rho)
P   <- matrix(c(1, rho1,
                rho1, 1), nrow = 2, byrow = TRUE)   # matriz das autocorrelacoes
rho <- c(rho1, rho2)
Phi <- solve(P) %*% rho                              # solve() = inversa/resolve sistema
cat("Pela matriz: phi12 =", round(Phi[1], 4), " phi22 =", round(Phi[2], 4), "\n")

# Conferencia 2: as 2 equacoes de Yule-Walker voltam aos dados?
cat("phi12 + rho1*phi22 =", round(phi12 + rho1*phi22, 4), "(deve dar", rho1, ")\n")
cat("rho1*phi12 + phi22 =", round(rho1*phi12 + phi22, 4), "(deve dar", rho2, ")\n")

# Conferencia 3 (extra, nao e do professor): PACF teorica pelo R
print(round(ARMAacf(ar = c(phi12, phi22), lag.max = 3, pacf = TRUE), 4))

# O AR(2) e estacionario? (triangulo dos slides)
cat("a1+a2<1:", phi12 + phi22 < 1, "| a2-a1<1:", phi22 - phi12 < 1,
    "| |a2|<1:", abs(phi22) < 1, "\n")
