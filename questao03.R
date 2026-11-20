# ==========================================================================
# 3o TRABALHO PRATICO - QUESTAO 3
# Considere o modelo AR(1):  Xt = 0,8*X(t-1) + Zt
#   a) O modelo e estacionario?
#   b) O modelo e inversivel?
# ==========================================================================

library(tseries)      # traz o kpss.test (o mesmo teste usado pelo professor)
graphics.off()        # fecha graficos antigos que possam ter ficado abertos
set.seed(13679821) # semente do professor (Aulas/09/ProcEstocasticos.r): trava o sorteio
                      # e os resultados saem IGUAIS toda vez que o script for rodado
n    <- 1000          # numero de valores da serie
alfa <- 0.8           # o numero do modelo: Xt = alfa*X(t-1) + Zt

# --------------------------------------------------------------------------
# (a) ESTACIONARIEDADE
# --------------------------------------------------------------------------
# Gera a serie como o professor faz no 01-AR.R:
# (aqui guardamos os sorteios Z no vetor 'zt' para usar no item (b))
zt <- rnorm(n, 0, 1)
xt <- 0
for (i in 2:n) {
  xt[i] <- alfa * xt[i - 1] + zt[i]  # Xt = 0,8*X(t-1) + Zt
}
plot(xt, type = "l", col = "red", main = "AR(1): Xt = 0,8 X(t-1) + Zt")

# Teste KPSS: H0 = serie estacionaria. p >= 0,05 -> NAO rejeita (e estacionaria).
kpss.test(xt, null = "L")

# Teoria (slide): raiz de Phi(B) = 1 - 0,8*B = 0 deve ficar FORA do circulo unitario
raiz <- polyroot(c(1, -alfa))        # coeficientes do menor para o maior grau: 1 e -0,8
raiz                                 # 1,25 (+0i)
Mod(raiz)                            # modulo = 1,25 > 1  -> FORA do circulo -> estacionario
abs(alfa) < 1                        # o atalho do AR(1): |alfa| < 1 -> TRUE

# Media e variancia teoricas (slide): E[X] = 0 ; Var[X] = sigma^2/(1 - alfa^2)
mean(xt)                             # perto de 0
var(xt)                              # perto de 2,78 (muda um pouco a cada sorteio)
1 / (1 - alfa^2)                     # variancia teorica = 2,78 (com sigma^2 = 1)

# ACF teorica (slide e 01-AR.R): rho(k) = alfa^k, para k = 0, 1, 2, ...
acft <- 0
for (k in 0:20) {
  acft[k + 1] <- alfa^k
}
acft                                 # 1, 0,8, 0,64, 0,512, ... decai aos poucos
lag <- 0:20
plot(lag, acft, type = "h", col = "blue", main = "ACF teorica do AR(1)")

# ACF empirica (estimada da serie gerada), comparada com a teorica:
a <- acf(xt, lag.max = 20)           # lag.max = 20 para ter o mesmo tamanho de 'acft'
plot(a$lag, a$acf, type = "h", col = "red", xlab = "Lag", ylab = "ACF")
lines(a$lag, acft, type = "b", col = "blue")   # azul = teorica, vermelho = estimada
# As duas ficam parecidas: a covariancia so depende do lag, nao do tempo.

pacf(xt, lag.max = 20)               # PACF: uma barra so no lag 1 (assinatura do AR(1))

# Conferencia com arima.sim (como o professor faz):
ar1 <- arima.sim(n = 1000, list(ar = c(alfa)))
par(mfrow = c(1, 3))
plot.ts(ar1, main = "AR(1)")
acf(ar1, type = "correlation", lag.max = 20, main = "AR(1)")
pacf(ar1, lag.max = 20, main = "AR(1)")
par(mfrow = c(1, 1))

# --------------------------------------------------------------------------
# (b) INVERSIBILIDADE
# --------------------------------------------------------------------------
# A condicao de inversibilidade (slides) e sobre theta(B), a parte de MEDIA MOVEL.
# No AR(1) puro nao ha parte MA: theta(B) = 1, que nao tem raiz alguma que reprove.
# Alem disso, a propria equacao ja entrega o Z com so dois termos:
#   Zt = Xt - 0,8*X(t-1)
z_recuperado <- xt[2:n] - alfa * xt[1:(n - 1)]
max(abs(z_recuperado - zt[2:n]))     # praticamente 0: recuperamos o Z exatamente a partir dos X

# --------------------------------------------------------------------------
# RESUMO
# --------------------------------------------------------------------------
data.frame(
  pergunta = c("(a) Estacionario?", "(b) Inversivel?"),
  resposta = c("SIM", "SIM"),
  motivo   = c("raiz de 1 - 0,8B = 0 e 1,25 > 1 (fora do circulo), equivale a |0,8| < 1",
               "todo AR e inversivel: Zt = Xt - 0,8*X(t-1), expressao finita (theta(B) = 1)")
)
