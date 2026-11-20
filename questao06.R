# ==========================================================================
# 3o TRABALHO PRATICO - QUESTAO 6
# Serie temporal:  Xt = (3/2)*X(t-1) - (1/2)*X(t-2) + Zt
#   a) Mostrar que a serie e NAO estacionaria.
#   b) Mostrar que Yt = diferenca de Xt (Yt = Xt - X(t-1)) e estacionario.
#   c) Simular no R uma serie de 1000 valores para Xt e, a partir dela,
#      produzir uma serie de 998 valores de Yt.
#   d) Plotar o correlograma de Yt e comentar o resultado.
#
# BASEADO NAS AULAS DO PROFESSOR:
#   Aulas/10/01-AR.R e Aulas-11/00-AR.R : gerar o AR(2) com o 'for' (xt <- 1:2 ...),
#       kpss.test, acf, pacf, arima.sim
#   Slides "Processos Autoregressivos": raizes de Phi(B) = 0 fora do circulo unitario
#       e as condicoes do triangulo do AR(2)
#   Slides Aulas-11 (ARMA-ARIMA, "Modelos Nao Estacionarios"):
#       Wt = (1 - B)^d Xt ; "se a serie e nao estacionaria o decaimento do correlograma
#       e muito lento" ; d = 1 e o valor mais usado
#   Aulas/09/ProcEstocasticos.r : diff() e acf() (passeio aleatorio)
# ==========================================================================

library(tseries)      # traz o kpss.test e o adf.test (testes de estacionariedade)
graphics.off()        # fecha graficos antigos que possam ter ficado abertos
set.seed(13679821)    # semente do professor (Aulas/09/ProcEstocasticos.r): trava o sorteio
                      # e os resultados saem IGUAIS toda vez que o script for rodado
n <- 1000             # numero de valores da serie Xt (o enunciado pede 1000)

# --------------------------------------------------------------------------
# (a) Xt e NAO estacionaria?  Teoria: raizes do polinomio Phi(B)
# --------------------------------------------------------------------------
# Xt = 1,5*X(t-1) - 0,5*X(t-2) + Zt     ->  alfa1 = 1,5  e  alfa2 = -0,5
# Phi(B) = 1 - 1,5*B + 0,5*B^2 = 0. As raizes precisam ficar FORA do circulo unitario.
raizes <- polyroot(c(1, -1.5, 0.5))  # coeficientes do menor para o maior grau: 1, -1,5, 0,5
raizes                               # 1 e 2
Mod(raizes)                          # modulos: 1 e 2. A raiz 1 esta SOBRE o circulo (raiz unitaria)

# Condicoes do triangulo para o AR(2) (as 3 precisam ser TRUE para ser estacionario):
a1 <- 1.5
a2 <- -0.5
a1 + a2 < 1                          # 1 < 1  -> FALSE  (reprovou, deu exatamente 1)
a2 - a1 < 1                          # -2 < 1 -> TRUE
(a2 > -1) & (a2 < 1)                 # -1 < -0,5 < 1 -> TRUE

# --------------------------------------------------------------------------
# (c) Simular 1000 valores de Xt (como o professor faz no 01-AR.R, exemplo do AR(2))
# --------------------------------------------------------------------------
xt <- c(0, 0)                        # dois valores iniciais (o AR(2) precisa de t-1 e t-2)
for (i in 3:n) {
  xt[i] <- 1.5 * xt[i - 1] - 0.5 * xt[i - 2] + rnorm(1, 0, 1)
}
length(xt)                           # 1000 valores
plot(xt, type = "l", col = "red", main = "Xt: nao estacionaria")
# A serie vaga e se afasta sem voltar para uma media fixa.

# Testes (os mesmos comentarios do professor no 00-AR.R):
kpss.test(xt, null = "L")            # H0: estacionaria.  p < 0,05 -> REJEITA -> nao estacionaria
adf.test(xt, alternative = "stationary")   # H0: raiz unitaria.  p > 0,05 -> NAO rejeita

acf(xt, lag.max = 20, main = "ACF de Xt")  # cai MUITO devagar = sinal de nao estacionariedade

# --------------------------------------------------------------------------
# (b) e (c) Yt = Xt - X(t-1)  (a primeira diferenca)
# --------------------------------------------------------------------------
# Teoria: Phi(B) = (1 - B)(1 - 0,5*B). O fator (1 - B) e a raiz unitaria, e e a diferenca:
#   Yt = (1 - B)*Xt = Xt - X(t-1)   ->   (1 - 0,5*B)*Yt = Zt   ->   Yt = 0,5*Y(t-1) + Zt
# Ou seja, Yt e um AR(1) com alfa = 0,5. A raiz de 1 - 0,5*B = 0 e B = 2 (> 1): ESTACIONARIO.
Mod(polyroot(c(1, -0.5)))            # 2 > 1  -> fora do circulo -> Yt e estacionario

# Calculo no R: diff(xt) tem 999 valores (Y2, ..., Y1000).
# O primeiro (Y2 = X2 - X1) usa os dois valores iniciais (0 e 0) que NOS escolhemos,
# entao descartamos ele e ficam os 998 valores que o enunciado pede.
yt <- diff(xt)[-1]
length(yt)                           # 998 valores

plot(yt, type = "l", col = "blue", main = "Yt = diferenca de Xt: estacionaria")
kpss.test(yt, null = "L")            # agora p >= 0,05: NAO rejeita -> estacionaria
adf.test(yt, alternative = "stationary")   # agora p < 0,05: rejeita a raiz unitaria -> estacionaria

# --------------------------------------------------------------------------
# (d) Correlograma de Yt (ACF e PACF)
# --------------------------------------------------------------------------
par(mfrow = c(1, 3))
plot(yt, type = "l", col = "blue", main = "Yt")
acf(yt, type = "correlation", lag.max = 20, main = "Correlograma (ACF) de Yt")
pacf(yt, lag.max = 20, main = "PACF de Yt")
par(mfrow = c(1, 1))

# ACF teorica de um AR(1) com alfa = 0,5 (como o professor faz no 01-AR.R): rho(k) = 0,5^k
acft <- 0
for (k in 0:20) {
  acft[k + 1] <- 0.5^k
}
acft                                 # 1, 0,5, 0,25, 0,125, ... decai rapido
lag <- 0:20
a <- acf(yt, lag.max = 20, plot = FALSE)
plot(a$lag, a$acf, type = "h", col = "red", xlab = "Lag", ylab = "ACF",
     main = "ACF de Yt (vermelho) x ACF teorica do AR(1) com alfa = 0,5 (azul)")
lines(a$lag, acft, type = "b", col = "blue")
round(a$acf[2:6], 3)                 # estimada: lags 1 a 5 (rho(1) perto de 0,5)
round(acft[2:6], 3)                  # teorica:  0,5  0,25  0,125  0,062  0,031

# Conferindo que Yt e mesmo um AR(1) com alfa perto de 0,5 (funcao arima, como no curso):
arima(yt, order = c(1, 0, 0), include.mean = FALSE)

# Conferencia com arima.sim (como o professor faz): ARIMA(1,1,0) com ar = 0,5
# (nota: o Xt do enunciado e exatamente um ARIMA(1,1,0))
x_sim <- arima.sim(n = 1000, list(order = c(1, 1, 0), ar = c(0.5)))
par(mfrow = c(1, 3))
plot.ts(x_sim, main = "ARIMA(1,1,0)")
acf(x_sim, lag.max = 20, main = "ACF de Xt")
acf(diff(x_sim), lag.max = 20, main = "ACF da diferenca")
par(mfrow = c(1, 1))

# --------------------------------------------------------------------------
# RESUMO
# --------------------------------------------------------------------------
data.frame(
  item = c("(a) Xt estacionaria?", "(b) Yt estacionaria?", "(d) correlograma de Yt"),
  resposta = c("NAO", "SIM", "decai rapido (AR(1), alfa = 0,5)"),
  motivo = c("raiz B = 1 sobre o circulo unitario; alfa1 + alfa2 = 1 (nao < 1)",
             "Yt = 0,5*Y(t-1) + Zt: raiz B = 2 > 1",
             "rho(k) = 0,5^k: 0,5; 0,25; 0,125...; PACF corta no lag 1")
)

# FRASE TECNICA: uma raiz unitaria (|B| = 1) torna a serie nao estacionaria e faz o
# correlograma decair muito devagar; uma diferenciacao (d = 1) elimina essa raiz e o
# Xt vira um ARIMA(p,1,q), cuja diferenca e um ARMA estacionario.
