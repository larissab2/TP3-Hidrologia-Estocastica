# ==========================================================================
# 3o TRABALHO PRATICO - QUESTAO 4
# Modelo: X(t) = 0,1*X(t-1) + 0,9*X(t-2) + Z(t)
#   a) Verificar as condicoes de estacionariedade
#   b) Se nao for estacionario, como torna-lo estacionario?
# ==========================================================================

library(tseries)      # traz o kpss.test e o adf.test (usados pelo professor)
graphics.off()        # fecha graficos antigos que possam ter ficado abertos
set.seed(13679821) # semente do professor (Aulas/09/ProcEstocasticos.r): trava o sorteio
                      # e os resultados saem IGUAIS toda vez que o script for rodado

# --------------------------------------------------------------------------
# (a) ESTACIONARIEDADE - parte teorica (as contas do trabalho)
# --------------------------------------------------------------------------
# Equacao caracteristica do AR(2):  1 - 0,1*B - 0,9*B^2 = 0
# polyroot() pede os numeros do termo sem B ate o B^2: 1, -0,1, -0,9
raizes <- polyroot(c(1, -0.1, -0.9))
raizes                # as raizes: 1 e -1,11
Mod(raizes)           # "tamanho" de cada raiz. Para ser estacionario TODAS > 1

# Condicoes do triangulo do AR(2) (as 3 precisam ser verdadeiras)
a1 <- 0.1
a2 <- 0.9
a1 + a2 < 1           # FALSE -> reprovou (deu exatamente 1)
a2 - a1 < 1           # TRUE
abs(a2) < 1           # TRUE

# --------------------------------------------------------------------------
# (a) ESTACIONARIEDADE - parte pratica: gerar a serie e testar (estilo do professor)
# --------------------------------------------------------------------------
# Gera X(t) = 0,1*X(t-1) + 0,9*X(t-2) + Z(t), com Z ~ N(0,1), 1000 valores
xt <- c(0, 0)                      # dois valores iniciais (o modelo precisa de t-1 e t-2)
for (i in 3:1000) {
  xt[i] <- 0.1*xt[i-1] + 0.9*xt[i-2] + rnorm(1, 0, 1)
}
plot(xt, type = "l", col = "red", main = "Serie X(t)")
# A serie "passeia" e nao volta para uma media fixa.

# Teste KPSS: H0 = a serie e estacionaria. Se o p-valor < 0,05, REJEITA H0.
kpss.test(xt, null = "L")
# Teste ADF: H0 = a serie TEM raiz unitaria. Se o p-valor > 0,05, NAO rejeita
# (ou seja, nao da para dizer que e estacionaria).
adf.test(xt, alternative = "stationary")

# ACF e PACF de X
acf(xt)                            # cai bem devagar = sinal de raiz unitaria
pacf(xt)

# --------------------------------------------------------------------------
# (b) COMO TORNAR ESTACIONARIO: diferenciar uma vez
# --------------------------------------------------------------------------
# Fatorando: 1 - 0,1B - 0,9B^2 = (1 - B)(1 + 0,9B)
# (1 - B) e a raiz unitaria -> a diferenca Y(t) = X(t) - X(t-1) elimina ela.
# Sobra (1 + 0,9B)*Y(t) = Z(t), ou seja, Y(t) = -0,9*Y(t-1) + Z(t): um AR(1)
Mod(polyroot(c(1, 0.9)))           # 1,11 > 1 -> Y e estacionario

yt <- diff(xt)                     # Y(t) = X(t) - X(t-1)
plot(yt, type = "l", col = "blue", main = "Serie Y(t) = diff(X)")

kpss.test(yt, null = "L")          # agora o p-valor deve ficar > 0,05 (nao rejeita)
adf.test(yt, alternative = "stationary")   # agora o p-valor deve ficar < 0,05

acf(yt)                            # alterna de sinal e some rapido (AR(1) com alfa = -0,9)
pacf(yt)                           # corta no lag 1 (assinatura do AR(1))

# --------------------------------------------------------------------------
# Conferencia com arima.sim (como o professor faz): ARIMA(1,1,0) com ar = -0,9
# --------------------------------------------------------------------------
arima110 <- arima.sim(n = 1000, list(order = c(1, 1, 0), ar = c(-0.9)))
par(mfrow = c(1, 3))
plot.ts(arima110, main = "ARIMA(1,1,0)")
acf(arima110, type = "correlation", lag.max = 20, main = "ARIMA(1,1,0)")
pacf(arima110, lag.max = 20, main = "ARIMA(1,1,0)")
par(mfrow = c(1, 1))

# Os 4 graficos principais juntos (uma figura so, boa para colocar no trabalho)
par(mfrow = c(2, 2))
plot(xt, type = "l", col = "red",  main = "Serie X(t)")
plot(yt, type = "l", col = "blue", main = "Serie Y(t) = diff(X)")
acf(xt, main = "ACF de X(t)")
acf(yt, main = "ACF de Y(t)")
par(mfrow = c(1, 1))

