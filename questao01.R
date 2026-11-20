# ==========================================================================
# 3o TRABALHO PRATICO - QUESTAO 1
# Os processos abaixo sao estacionarios? Justifique.   Zt ~ iid N(0, sigma^2)
#   (a) Xt = Zt - Z(t-1)
#   (b) Xt = X(t-1) + Zt
#   (c) Xt - X(t-1) + 0,5*X(t-2) = Zt
# ==========================================================================

library(tseries)      # traz o kpss.test (o mesmo teste usado pelo professor)
graphics.off()        # fecha graficos antigos que possam ter ficado abertos
set.seed(13679821) # semente do professor (Aulas/09/ProcEstocasticos.r): trava o sorteio
                      # e os resultados saem IGUAIS toda vez que o script for rodado
n <- 1000             # numero de valores de cada serie (como no 01-MA.R)

# --------------------------------------------------------------------------
# (a) Xt = Zt - Z(t-1)   ->   MA(1) com beta1 = -1
# --------------------------------------------------------------------------
# Gera a serie como o professor faz no 01-MA.R:
zt <- rnorm(1, 0, 1)                 # primeiro sorteio
xa <- 0
for (i in 2:(n + 1)) {
  zt[i] <- rnorm(1, 0, 1)            # sorteio novo (o "Z" de hoje)
  xa[i] <- zt[i] - 1 * zt[i - 1]     # Xt = Zt + beta1*Z(t-1), com beta1 = -1
}
xa <- xa[2:(n + 1)]                  # descarta o primeiro valor (so um zero inicial)

mean(xa)                             # media: deve ficar perto de 0 (constante)
var(xa)                              # variancia: perto de 2 (= 2*sigma^2); muda um pouco a cada sorteio

# Teste KPSS: H0 = serie estacionaria. Se p >= 0,05, NAO rejeita (e estacionaria).
kpss.test(xa, null = "L")

# Valores teoricos (conta feita a mao): gama(0) = 2, gama(1) = -1, rho(1) = -0,5
# Pela formula do slide do MA(1): rho(1) = beta1/(1 + beta1^2) = -1/(1+1) = -0,5
beta1 <- -1
beta1 / (1 + beta1^2)                # rho(1) teorico = -0,5
acf(xa, lag.max = 20, plot = FALSE)$acf[2]   # rho(1) estimado dos dados (perto de -0,5)

par(mfrow = c(1, 3))
plot(xa, type = "l", col = "red", main = "(a) Xt = Zt - Z(t-1)")
acf(xa, type = "correlation", lag.max = 20, main = "(a) ACF")
pacf(xa, lag.max = 20, main = "(a) PACF")
par(mfrow = c(1, 1))
# A serie oscila em torno de 0 sem se afastar, e a ACF so tem barra no lag 1.

# Conferencia com arima.sim (como o professor faz). Aqui ma = c(-1).
ma_a <- arima.sim(n = 1000, list(ma = c(-1)))
var(ma_a)                            # tambem perto de 2

# --------------------------------------------------------------------------
# (b) Xt = X(t-1) + Zt   ->   PASSEIO ALEATORIO (aula 09)
# --------------------------------------------------------------------------
# Gera a serie como o professor faz no ProcEstocasticos.r:
xb <- 0
for (i in 2:n) {
  xb[i] <- xb[i - 1] + rnorm(1, mean = 0, sd = 1)   # Xt = X(t-1) + Zt
}
kpss.test(xb, null = "L")            # p-valor pequeno (< 0,05): REJEITA H0 -> NAO estacionaria

par(mfrow = c(2, 2))                 # 4 graficos juntos, como no ProcEstocasticos.r
plot(xb, type = "l", xlab = "Tempo", ylab = "Xt", main = "(b) Passeio aleatorio")
acf(xb, main = "ACF de Xt", xlab = "Lag", ylab = "autocorrelacoes")
zb <- diff(xb)                       # primeira diferenca: Xt - X(t-1) = Zt
plot(zb, type = "l", xlab = "Tempo", ylab = "1a diferenca", main = "1a diferenca de Xt")
acf(zb, main = "ACF da 1a diferenca", xlab = "Lag", ylab = "autocorrelacoes")
par(mfrow = c(1, 1))
# A serie vaga sem voltar a media e a ACF quase nao cai. A diferenca e ruido branco.

kpss.test(zb, null = "L")            # agora p >= 0,05: a diferenca e estacionaria

# Por que NAO e estacionario? Var(Xt) = t * sigma^2 CRESCE com o tempo (slide da aula 09).
# Conferimos sorteando 2000 "caminhos" de 100 passos e olhando a variancia em cada t:
caminhos <- 2000
tmax     <- 100
Z <- matrix(rnorm(caminhos * tmax), nrow = caminhos, ncol = tmax)  # cada linha = um caminho
X <- t(apply(Z, 1, cumsum))          # cumsum: Xt = Z1 + Z2 + ... + Zt (comecando em X0 = 0)
tempos <- c(1, 2, 3, 10, 50, 100)
tabela <- data.frame(
  tempo        = tempos,
  var_simulada = apply(X[, tempos], 2, var),   # variancia de verdade, medida nos caminhos
  var_teorica  = tempos * 1                    # t * sigma^2, com sigma^2 = 1
)
tabela                               # a variancia cresce junto com t: NAO estacionario

# --------------------------------------------------------------------------
# (c) Xt - X(t-1) + 0,5*X(t-2) = Zt   ->   Xt = 1*X(t-1) - 0,5*X(t-2) + Zt   (AR(2))
#     alfa1 = 1  e  alfa2 = -0,5
# --------------------------------------------------------------------------
# Gera a serie como o professor faz no 01-AR.R (exemplo do AR(2)):
xc <- c(0, 0)                        # dois valores iniciais
for (i in 3:n) {
  xc[i] <- 1 * xc[i - 1] - 0.5 * xc[i - 2] + rnorm(1, 0, 1)
}
kpss.test(xc, null = "L")            # p >= 0,05: nao rejeita -> estacionaria

# Teoria (slide "Estacionariedade de um processo AR(p)"):
# polinomio Phi(B) = 1 - alfa1*B - alfa2*B^2 = 1 - B + 0,5*B^2 ; precisa ter raizes FORA do circulo
raizes <- polyroot(c(1, -1, 0.5))    # coeficientes do menor para o maior grau: 1, -1, 0,5
raizes                               # 1+1i e 1-1i (numeros complexos)
Mod(raizes)                          # modulo = 1,414 > 1  -> raizes FORA do circulo unitario

# Condicoes do triangulo para o AR(2) (as 3 precisam ser TRUE):
a1 <- 1
a2 <- -0.5
a1 + a2 < 1                          # 0,5 < 1   -> TRUE
a2 - a1 < 1                          # -1,5 < 1  -> TRUE
(a2 > -1) & (a2 < 1)                 # -1 < -0,5 < 1 -> TRUE

par(mfrow = c(1, 3))
plot(xc, type = "l", col = "blue", main = "(c) AR(2)")
acf(xc, lag.max = 20, main = "(c) ACF")
pacf(xc, lag.max = 20, main = "(c) PACF")
par(mfrow = c(1, 1))
# Conferencia com arima.sim (como o professor faz)
ar2 <- arima.sim(n = 1000, list(ar = c(1, -0.5)))
par(mfrow = c(1, 3))
plot.ts(ar2, main = "AR(2)")
acf(ar2, type = "correlation", lag.max = 20, main = "AR(2)")
pacf(ar2, lag.max = 20, main = "AR(2)")
par(mfrow = c(1, 1))

# --------------------------------------------------------------------------
# RESUMO
# --------------------------------------------------------------------------
data.frame(
  processo      = c("(a) Xt = Zt - Z(t-1)", "(b) Xt = X(t-1) + Zt", "(c) AR(2) com alfa1=1 e alfa2=-0,5"),
  estacionario  = c("SIM", "NAO", "SIM"),
  motivo        = c("media 0, variancia 2*sigma^2 e covariancia so depende do lag",
                    "Var(Xt) = t*sigma^2 cresce com t (raiz unitaria)",
                    "raizes 1+-1i com modulo 1,414 > 1 (fora do circulo)")
)
