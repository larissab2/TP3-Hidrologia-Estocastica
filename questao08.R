# =====================================================================
# 3o TRABALHO PRATICO - HIDROLOGIA ESTOCASTICA (EHR806)
# QUESTAO 8 - AR(2) estacionario: Xt = a1 Xt-1 + a2 Xt-2 + Zt
#             Sabendo a1 = 0,7 e rho2 = 0,6, da para achar a2?
# =====================================================================
graphics.off()
alfa1 <- 0.7
rho2  <- 0.6

# Montando a equacao: a1^2/(1-a2) + a2 = rho2.
# Multiplicando por (1-a2) e arrumando fica uma equacao de 2o grau em a2:
#     a2^2 - (1 + rho2)*a2 + (rho2 - alfa1^2) = 0
A <- 1
B <- -(1 + rho2)          # = -1,6
C <- rho2 - alfa1^2       # =  0,11
delta <- B^2 - 4*A*C      # discriminante (Bhaskara)
a2 <- c((-B - sqrt(delta))/(2*A), (-B + sqrt(delta))/(2*A))   # as 2 raizes
cat("Duas raizes da equacao:", round(a2, 4), "\n")

# Para cada raiz, testo as 3 condicoes de estacionariedade do AR(2) (triangulo)
for (x in a2) {
  c1 <- (alfa1 + x) < 1
  c2 <- (x - alfa1) < 1
  c3 <- abs(x) < 1
  cat("alfa2 =", round(x, 4), "| a1+a2<1:", c1, "| a2-a1<1:", c2,
      "| |a2|<1:", c3, "=> estacionario:", c1 & c2 & c3, "\n")
}

# Fico com a raiz que passa nos 3 testes e confiro de volta
a2.ok <- a2[abs(a2) < 1]
rho1  <- alfa1/(1 - a2.ok)
cat("alfa2 =", round(a2.ok, 4), " rho1 =", round(rho1, 4), "\n")
cat("Conferencia: rho2 = a1*rho1 + a2 =", round(alfa1*rho1 + a2.ok, 4), "\n")

# Conferencia extra: ACF teorica do AR(2) encontrado, no estilo do professor
# (Aulas/10/01-AR.R calcula a ACF teorica com um for e plota com type="h").
# Formula dos slides do AR(2): rho(k) = a1*rho(k-1) + a2*rho(k-2)
acft <- numeric(21)
acft[1] <- 1                    # rho(0) = 1
acft[2] <- rho1                 # rho(1) que achamos acima
for (k in 2:20) {
  acft[k+1] <- alfa1*acft[k] + a2.ok*acft[k-1]
}
print(round(acft[1:4], 4))      # deve dar 1 ; 0,7543 ; 0,6 ; 0,4743
lag <- 0:20
plot(lag, acft, type = "h", col = "blue", main = "ACF teorica do AR(2) encontrado")
abline(h = 0)

# Raizes do polinomio 1 - a1 B - a2 B^2: devem ter modulo > 1 (estacionario).
# (polyroot nao aparece nos scripts do professor; e so uma checagem extra)
print(round(Mod(polyroot(c(1, -alfa1, -a2.ok))), 4))
