# =====================================================================
# 3o TRABALHO PRATICO - HIDROLOGIA ESTOCASTICA (EHR806)
# QUESTAO 9 - Determinar a ACF de  Xt = 14 + Zt + 0,4 Zt-1 - 0,2 Zt-2
# Baseado em: Aulas/10/00-Processos de Media Movel.pdf
#     (rho(k) do MA(q); somar uma constante mu NAO muda a ACF)
#             Aulas/10/00-MA.R (serie gerada com for + acf)
# =====================================================================
library(tseries)
graphics.off()
set.seed(13679821)       # semente do professor (Aulas/09/ProcEstocasticos.r)

# ---------------------------------------------------------------------
# PARTE 1 - ACF TEORICA (e a resposta da questao)
# O 14 e so uma constante (mu): NAO entra na ACF. Fica um MA(2):
#   Xt = Zt + beta1*Zt-1 + beta2*Zt-2  com beta0 = 1
# ---------------------------------------------------------------------
beta <- c(1, 0.4, -0.2)            # beta0, beta1, beta2
q <- 2
soma <- sum(beta^2)                # 1 + 0,16 + 0,04 = 1,2  (e o "rho(0)" sem dividir)
cat("soma dos beta^2 =", soma, "\n")

rho <- numeric(21)                 # lags 0 a 20
rho[1] <- 1                        # rho(0) = 1
for (k in 1:20) {
  if (k <= q) {
    rho[k+1] <- sum(beta[1:(q-k+1)] * beta[(1+k):(q+1)]) / soma
  } else {
    rho[k+1] <- 0                  # depois do lag q, zero (corta)
  }
}
cat("rho(1) =", round(rho[2], 4), "  (= 0,32/1,2)\n")
cat("rho(2) =", round(rho[3], 4), "  (= -0,2/1,2)\n")
cat("rho(k>=3) =", rho[4:6], "\n")

lag <- 0:20
plot(lag, rho, type = "h", col = "blue", ylim = c(-1, 1),
     main = "ACF teorica de Xt = 14 + Zt + 0,4Zt-1 - 0,2Zt-2")
abline(h = 0)

# ---------------------------------------------------------------------
# PARTE 2 - CONFERENCIA POR SIMULACAO (como no 00-MA.R do professor)
# Se a teoria estiver certa, a ACF da serie simulada fica parecida.
# ---------------------------------------------------------------------
n <- 5000
zt <- rnorm(n + 2, 0, 1)           # os choques Z (ruido branco)
xt <- numeric(n)
for (i in 1:n) {
  xt[i] <- 14 + zt[i+2] + 0.4*zt[i+1] - 0.2*zt[i]   # a formula do enunciado
}
cat("media da serie simulada =", round(mean(xt), 3), "(deve ser perto de 14)\n")
a <- acf(xt, lag.max = 20, plot = FALSE)
cat("ACF simulada lag 1 a 3:", round(a$acf[2:4], 3), "\n")
par(mfrow = c(1, 2))
plot(xt[1:200], type = "l", col = "red", main = "Serie simulada (200 primeiros)")
acf(xt, type = "correlation", lag.max = 20, main = "ACF da serie simulada")
par(mfrow = c(1, 1))

# FRASE TECNICA: a constante 14 apenas desloca a media (E[Xt] = 14) e nao afeta
# a ACF; o processo e um MA(2) com rho(1) = 0,32/1,2 = 0,267, rho(2) = -0,2/1,2 = -0,167
# e rho(k) = 0 para k > 2 (a ACF corta no lag q = 2).
