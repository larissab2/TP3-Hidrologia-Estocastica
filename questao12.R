# =====================================================================
# 3o TRABALHO PRATICO - HIDROLOGIA ESTOCASTICA (EHR806)
# QUESTAO 12 - n = 100; rho1=0,75, rho2=0,09, rho3=0,04, rho4=0,02, rho5=0,015
#   a) Calcular phi11, phi22, phi33 pelas equacoes de Yule-Walker
# =====================================================================
graphics.off()
n   <- 100
rho <- c(0.75, 0.09, 0.04, 0.02, 0.015)    # rho1 a rho5 (amostrais)

# ---------------------------------------------------------------------
# 1) Pela forma matricial do slide (eq. 05):  Phi_k = inversa(P_k) %*% rho_k
# ---------------------------------------------------------------------
phi.kk <- numeric(3)
for (k in 1:3) {
  P <- matrix(1, k, k)                      # matriz das autocorrelacoes
  for (i in 1:k) for (j in 1:k) {
    if (i != j) P[i, j] <- rho[abs(i - j)]  # rho(|i-j|) fora da diagonal
  }
  Phi <- solve(P) %*% rho[1:k]              # resolve o sistema
  phi.kk[k] <- Phi[k]                       # o ULTIMO coeficiente e o phi_kk
}
cat("phi11, phi22, phi33 (matriz) =", round(phi.kk, 4), "\n")

# ---------------------------------------------------------------------
# 2) Conferencia pelas relacoes de Durbin (eqs. 06, 07, 08, 09 dos slides)
# ---------------------------------------------------------------------
f11 <- rho[1]                                         # eq. 06
f22 <- (rho[2] - rho[1]^2)/(1 - rho[1]^2)             # eq. 07
f12 <- f11 - f22*f11                                  # eq. 09 (k=2, i=1)
num <- rho[3] - (f12*rho[2] + f22*rho[1])             # eq. 08 (k=3): numerador
den <- 1      - (f12*rho[1] + f22*rho[2])             #               denominador
f33 <- num/den
cat("phi11, phi22, phi33 (Durbin) =", round(c(f11, f22, f33), 4), "\n")
cat("  (phi12 =", round(f12, 4), "| numerador =", round(num, 4),
    "| denominador =", round(den, 4), ")\n")

# ---------------------------------------------------------------------
# 3) ESTES NUMEROS FAZEM SENTIDO? (um coeficiente de correlacao nao passa de 1)
# ---------------------------------------------------------------------
cat("Algum |phi_kk| > 1 ?", any(abs(phi.kk) > 1), "\n")
cat("Dado rho1 =", rho[1], "o rho2 teria que ser >=", 2*rho[1]^2 - 1,
    "(e o enunciado traz", rho[2], ")\n")
P3 <- matrix(c(1, rho[1], rho[2],
               rho[1], 1, rho[1],
               rho[2], rho[1], 1), 3, 3, byrow = TRUE)
cat("Determinante da matriz de autocorrelacoes (deveria ser > 0):", round(det(P3), 4), "\n")

# Faixa de confianca 95% usada no R (acf/pacf): +- 1,96/sqrt(n)
cat("Faixa +-", round(1.96/sqrt(n), 3), "\n")
