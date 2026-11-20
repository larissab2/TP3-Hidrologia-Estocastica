# ==========================================================================
# 3o TRABALHO PRATICO - QUESTAO 2
# Considere o modelo MA(1):  Xt = Zt - 0,9*Z(t-1)
#   a) O modelo e estacionario?
#   b) O modelo e inversivel?
# ==========================================================================

library(tseries)      # traz o kpss.test (o mesmo teste usado pelo professor)
graphics.off()        # fecha graficos antigos que possam ter ficado abertos
set.seed(13679821) # semente do professor (Aulas/09/ProcEstocasticos.r): trava o sorteio
                      # e os resultados saem IGUAIS toda vez que o script for rodado
n     <- 1000         # numero de valores da serie (como no 01-MA.R)
beta1 <- -0.9         # o numero do modelo: Xt = Zt + beta1*Z(t-1) = Zt - 0,9*Z(t-1)

# --------------------------------------------------------------------------
# (a) ESTACIONARIEDADE
# --------------------------------------------------------------------------
# Gera a serie como o professor faz no 01-MA.R:
zt <- rnorm(1, 0, 1)                 # primeiro sorteio
xt <- 0
for (i in 2:(n + 1)) {
  zt[i] <- rnorm(1, 0, 1)            # sorteio novo (o "Z" de hoje)
  xt[i] <- zt[i] + beta1 * zt[i - 1] # Xt = Zt - 0,9*Z(t-1)
}
xt <- xt[2:(n + 1)]                  # descarta o primeiro valor (so um zero inicial)
z_verdadeiro <- zt[2:(n + 1)]        # guardamos os sorteios reais, para usar no item (b)

# Valores teoricos (as contas do trabalho), com sigma^2 = 1:
#   media = 0 ; variancia = (1 + beta1^2) = 1,81 ; covariancia lag 1 = beta1 = -0,9
mean(xt)                             # perto de 0
var(xt)                              # perto de 1,81
1 + beta1^2                          # variancia teorica = 1,81
beta1 / (1 + beta1^2)                # rho(1) teorico = -0,497 (formula do slide)
acf(xt, lag.max = 20, plot = FALSE)$acf[2]  # rho(1) estimado dos dados (perto de -0,497)

# Teste KPSS: H0 = serie estacionaria. p >= 0,05 -> NAO rejeita (e estacionaria).
kpss.test(xt, null = "L")

par(mfrow = c(1, 3))
plot(xt, type = "l", col = "red", main = "MA(1): Xt = Zt - 0,9 Z(t-1)")
acf(xt, type = "correlation", lag.max = 20, main = "ACF")
pacf(xt, lag.max = 20, main = "PACF")
par(mfrow = c(1, 1))
# A ACF tem so uma barra (lag 1, negativa) e corta: assinatura do MA(1).

# Conferencia com arima.sim (como o professor faz):
ma1 <- arima.sim(n = 1000, list(ma = c(-0.9)))
par(mfrow = c(1, 3))
plot.ts(ma1, main = "MA(1)")
acf(ma1, type = "correlation", lag.max = 20, main = "MA(1)")
pacf(ma1, lag.max = 20, main = "MA(1)")
par(mfrow = c(1, 1))

# --------------------------------------------------------------------------
# (b) INVERSIBILIDADE
# --------------------------------------------------------------------------
# Polinomio do MA: theta(B) = 1 + beta1*B = 1 - 0,9*B. Precisa ter a raiz FORA do circulo.
raiz <- polyroot(c(1, beta1))        # coeficientes do menor para o maior grau: 1 e -0,9
raiz                                 # 1,11 (+0i)
Mod(raiz)                            # modulo = 1,11 > 1  -> FORA do circulo -> INVERSIVEL
abs(beta1) < 1                       # o atalho do MA(1): |beta1| < 1  -> TRUE

# O que "inversivel" quer dizer: dar para recuperar o Z a partir dos X.
# Z(t) = X(t) + 0,9*X(t-1) + 0,81*X(t-2) + ...  (pesos 0,9^j, que diminuem)
pesos <- (-beta1)^(0:10)
pesos                                # 1, 0,9, 0,81, 0,729, ... (vao diminuindo)

# Teste pratico: reconstruir os Z usando SO os X (somando os 50 primeiros termos)
z_estimado <- stats::filter(xt, filter = (-beta1)^(0:50), method = "convolution", sides = 1)
z_estimado <- as.numeric(z_estimado)
ok <- 51:n                           # ignora o comeco (faltam valores do passado)
cor(z_estimado[ok], z_verdadeiro[ok])   # perto de 1: recuperamos o Z a partir dos X
plot(z_verdadeiro[ok], z_estimado[ok], xlab = "Z verdadeiro (sorteado)",
     ylab = "Z recuperado dos X", main = "Z recuperado x Z verdadeiro")
abline(0, 1, col = "red")            # pontos em cima da reta = recuperou direitinho

# Por que exigimos |beta1| < 1? Porque beta1 = -0,9 e beta1 = -1/0,9 = -1,11 geram
# EXATAMENTE a mesma ACF, e o R nao saberia qual escolher:
beta_b <- -1 / 0.9
beta1 / (1 + beta1^2)                # rho(1) com beta1 = -0,9   = -0,497
beta_b / (1 + beta_b^2)              # rho(1) com beta1 = -1,11  = -0,497 (igual!)
Mod(polyroot(c(1, beta_b)))          # 0,9 < 1 -> este NAO e inversivel (so um dos dois e)

# --------------------------------------------------------------------------
# RESUMO
# --------------------------------------------------------------------------
data.frame(
  pergunta = c("(a) Estacionario?", "(b) Inversivel?"),
  resposta = c("SIM", "SIM"),
  motivo   = c("todo MA e estacionario: media 0, variancia 1,81*sigma^2, covariancia so depende do lag",
               "raiz de 1 - 0,9B = 0 e 1,11 > 1 (fora do circulo unitario), equivale a |0,9| < 1")
)

