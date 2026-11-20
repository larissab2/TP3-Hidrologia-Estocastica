# =====================================================================
# 3o TRABALHO PRATICO - HIDROLOGIA ESTOCASTICA (EHR806)
# QUESTAO 13 - n = 200. Dadas a FAC e a FACP estimadas, identificar um
#              possivel modelo gerador dos dados.
# Baseado em: Aulas-11/03-ARMA-ARIMA.pdf (tabela de identificacao ACF x PACF:
#             AR(p), MA(q), ARMA(p,q); "ordem do AR = numero de PACF
#             significativas"; "MA = maior lag significativo da ACF")
#             Aulas-11/02-ARMA.R (arima.sim + plot.ts + acf + pacf)
# =====================================================================
library(tseries)
graphics.off()
set.seed(13679821)       # semente do professor (Aulas/09/ProcEstocasticos.r)

n    <- 200
k    <- 1:5
rho  <- c(-0.909, 0.727, -0.54, 0.454, -0.364)    # FAC estimada
phi  <- c(-0.909, -0.317, -0.225, -0.136, -0.410) # FACP estimada
lim  <- 1.96/sqrt(n)                               # banda de 95% (padrao do R)
cat("Banda de 95%: +-", round(lim, 3), "\n")

# --- 1) DESENHAR a FAC e a FACP dadas (barrinhas + banda) ----------------
par(mfrow = c(1, 2))
plot(k, rho, type = "h", lwd = 4, col = "blue", ylim = c(-1, 1),
     main = "FAC estimada", xlab = "lag k", ylab = "rho(k)")
abline(h = 0); abline(h = c(-lim, lim), col = "red", lty = 2)
plot(k, phi, type = "h", lwd = 4, col = "darkgreen", ylim = c(-1, 1),
     main = "FACP estimada", xlab = "lag k", ylab = "phi(kk)")
abline(h = 0); abline(h = c(-lim, lim), col = "red", lty = 2)
par(mfrow = c(1, 1))

# --- 2) O QUE A FAC DIZ: decai ou corta? --------------------------------
# Dividindo cada rho pelo anterior vemos se "cada passo multiplica pelo mesmo numero"
razao <- rho[2:5]/rho[1:4]
cat("rho(k)/rho(k-1) =", round(razao, 3), " media =", round(mean(razao), 3), "\n")
# razao quase constante e NEGATIVA (~ -0,8): decai exponencialmente ALTERNANDO sinal
# => nao corta => nao e MA puro.

# --- 3) O QUE A FACP DIZ: corta no lag 1 ou decai? ----------------------
cat("FACP fora da banda (|phi| >", round(lim, 3), "):", abs(phi) > lim, "\n")
# lags 2 e 3 tambem passam da banda e vao diminuindo (-0,317; -0,225; -0,136)
# => nao e um corte limpo em p = 1 => tambem decai.
# Se fosse AR(1) com alfa = rho(1), a FAC seria alfa^k:
cat("AR(1) com alfa = -0,909 daria rho =", round((-0.909)^k, 3), "\n")
cat("  mas a FAC dada e                 ", rho, "\n")

# --- 4) ENTAO: FAC decai + FACP decai  ==>  ARMA(1,1) (tabela do professor) -
# alfa1 ~ razao media (~ -0,8). Procuro beta1 que reproduz rho(1) = -0,909:
# rho(1) = (1 + a*b)(a + b)/(1 + 2*a*b + b^2)   (formula classica Box-Jenkins)
alfa1 <- round(mean(razao), 2)
melhor <- c(NA, Inf)
for (b in seq(-0.99, 0.99, 0.01)) {
  r1 <- (1 + alfa1*b)*(alfa1 + b)/(1 + 2*alfa1*b + b^2)
  erro <- abs(r1 - rho[1])
  if (erro < melhor[2]) melhor <- c(b, erro)
}
beta1 <- melhor[1]
cat("Candidato: ARMA(1,1) com alfa1 =", alfa1, " e beta1 =", beta1, "\n")

# --- 5) CONFERIR: simular esse ARMA(1,1) e ver se a FAC/FACP se parecem ---
# (como no 02-ARMA.R do professor: arima.sim)
arma1 <- arima.sim(n = n, list(ar = c(alfa1), ma = c(beta1)))
par(mfrow = c(1, 3))
plot.ts(arma1, main = "ARMA(1,1) simulado")
acf(arma1, type = "correlation", lag.max = 20, main = "FAC simulada")
pacf(arma1, lag.max = 20, main = "FACP simulada")
par(mfrow = c(1, 1))
kpss.test(as.numeric(arma1), null = "L")   # H0: estacionaria

# --- 6) AVISO: a FACP da tabela e coerente com a FAC da tabela? ---------
# Pelas equacoes de Yule-Walker (Aula 10), a PACF decorre da FAC:
pac <- numeric(5)
for (kk in 1:5) {
  P <- diag(kk)
  for (i in 1:kk) for (j in 1:kk) if (i != j) P[i, j] <- rho[abs(i - j)]
  pac[kk] <- solve(P) %*% rho[1:kk] |> as.numeric() |> tail(1)
}
cat("FACP recalculada a partir da FAC:", round(pac, 3), "\n")
cat("FACP da tabela                  :", phi, "\n")
cat("-> nao batem (e passam de 1 nos lags finais): os dados do enunciado sao\n")
cat("   aproximados/inconsistentes; a identificacao e feita pelo PADRAO.\n")

# FRASE TECNICA: a FAC decai exponencialmente alternando de sinal (rho(k) ~ -0,8
# rho(k-1)) e a FACP tambem decai (nao ha corte nitido); pela tabela de
# identificacao de Box-Jenkins (ambas decaem) o modelo candidato e um ARMA(1,1)
# com alfa1 ~ -0,8 (alfa1 < 0 explica a alternancia de sinal).
