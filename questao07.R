# =====================================================================
# 3o TRABALHO PRATICO - HIDROLOGIA ESTOCASTICA (EHR806)
# QUESTAO 7 - ACF teorica: (a) AR(1) com alfa = 0,1; -0,75; 0,99
#                          (b) Medias moveis de ordem q
# =====================================================================
graphics.off()          # fecha graficos antigos (evita Plots vazio no RStudio)
lag.max <- 20           # vamos olhar do lag 0 ate o lag 20
k <- 0:lag.max          # k = "quantos passos para tras"

# ---------------------------------------------------------------------
# (a) AR(1): a formula do professor e  rho(k) = alfa^k
#     (cada passo para tras multiplica por alfa de novo)
# ---------------------------------------------------------------------
alfas <- c(0.1, -0.75, 0.99)
par(mfrow = c(1, 3))     # 3 graficos lado a lado
for (a in alfas) {
  rho <- a^k                                  # a formula do professor
  rho.R <- ARMAacf(ar = a, lag.max = lag.max) # o proprio R calcula (conferencia)
  cat("alfa =", a, "| maior diferenca formula x R =",
      format(max(abs(rho - rho.R)), digits = 3), "\n")
  print(round(rho[1:6], 4))                   # lags 0 a 5
  plot(k, rho, type = "h", lwd = 3, ylim = c(-1, 1), # type="h" = barrinhas
       xlab = "lag k", ylab = expression(rho(k)),
       main = paste("AR(1), alfa =", a))
  abline(h = 0)
}
par(mfrow = c(1, 1))     # volta ao normal (forma correta: mfrow = c(1,1))

# ---------------------------------------------------------------------
# (b) MA(q): rho(k) CORTA no lag q  (rho(k) = 0 para k > q)
#     Exemplo: MA(2)  Xt = Zt + 0,6 Zt-1 + 0,3 Zt-2   (beta0 = 1)
# ---------------------------------------------------------------------
beta <- c(1, 0.6, 0.3)                 # beta0, beta1, beta2
q <- length(beta) - 1                  # q = 2
rho.ma <- numeric(lag.max + 1); rho.ma[1] <- 1     # rho(0) = 1
for (kk in 1:lag.max) {
  if (kk <= q) {                       # formula do professor para k = 1..q
    rho.ma[kk + 1] <- sum(beta[1:(q - kk + 1)] * beta[(1 + kk):(q + 1)]) / sum(beta^2)
  } else {
    rho.ma[kk + 1] <- 0                # k > q  -> zero (corta)
  }
}
rho.ma.R <- ARMAacf(ma = beta[-1], lag.max = lag.max)
cat("MA(2): maior diferenca formula x R =",
    format(max(abs(rho.ma - rho.ma.R)), digits = 3), "\n")
print(round(rho.ma[1:6], 4))
plot(k, rho.ma, type = "h", lwd = 3, ylim = c(-1, 1),
     xlab = "lag k", ylab = expression(rho(k)),
     main = "MA(2): a ACF corta no lag q = 2")
abline(h = 0)
