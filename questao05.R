# ==========================================================================
# 3o TRABALHO PRATICO - QUESTAO 5
# AR(1): X(t) = alfa*X(t-1) + Z(t), escrito como MA infinito. Mostrar que
#   E[X(t)] = 0   e   Var[X(t)] = sigma^2 * (1 + alfa^2 + alfa^4 + ...) = sigma^2/(1 - alfa^2)
#
# Baseado nos slides do professor "Processos Autoregressivos" (Aulas/10):
#   substituicoes sucessivas -> X(t) = soma de alfa^(j-1) * Z(t-j+1), j = 1, 2, ...
#   e na forma com o operador B: X(t) = (1 + alfa*B + alfa^2*B^2 + ...) Z(t)
# Estilo do script do professor (Aulas/10/01-AR.R): for(i in ...) para gerar o
# AR e arima.sim. Aqui o script CONFERE numericamente a conta feita a mao.
# ==========================================================================

graphics.off()        # fecha graficos antigos que possam estar abertos
set.seed(13679821) # semente do professor (Aulas/09/ProcEstocasticos.r): trava o sorteio
                      # e os resultados saem IGUAIS toda vez que o script for rodado

alfa   <- 0.8         # coeficiente do AR(1) (pode trocar por outro, com |alfa| < 1)
sigma2 <- 1           # variancia do ruido Z (Z ~ N(0,1), como nas aulas)

# --------------------------------------------------------------------------
# 1) Os pesos do MA infinito: X(t) = Z(t) + alfa*Z(t-1) + alfa^2*Z(t-2) + ...
#    No slide: peso de Z(t-j+1) = alfa^(j-1), j = 1, 2, ...  (mesmos numeros que
#    alfa^j com j = 0, 1, 2, ...). Tambem sao os coeficientes de 1/(1 - alfa*B).
# --------------------------------------------------------------------------
j     <- 0:20
pesos <- alfa^j
pesos                           # 1, 0,8, 0,64, 0,512 ... (vao diminuindo)
plot(j, pesos, type = "h", col = "blue", xlab = "j (quantos passos atras)",
     ylab = "peso do Z(t-j)")

# A funcao ARMAtoMA do R calcula esses mesmos pesos (a partir de j = 1)
ARMAtoMA(ar = alfa, ma = 0, lag.max = 10)

# --------------------------------------------------------------------------
# 2) A variancia como soma: sigma^2 * (1 + alfa^2 + alfa^4 + ...)
# --------------------------------------------------------------------------
soma_ate_5   <- sigma2 * sum(alfa^(2 * (0:5)))     # so os 6 primeiros termos
soma_ate_50  <- sigma2 * sum(alfa^(2 * (0:50)))    # 51 termos
soma_ate_500 <- sigma2 * sum(alfa^(2 * (0:500)))   # quase infinitos
formula      <- sigma2 / (1 - alfa^2)              # o resultado pedido

c(ate_5_termos = soma_ate_5, ate_50_termos = soma_ate_50,
  ate_500_termos = soma_ate_500, formula = formula)
# A soma vai se aproximando de 2,78, que e o valor de sigma^2/(1 - alfa^2).

# --------------------------------------------------------------------------
# 3) Conferencia por simulacao: gerar um AR(1) com o 'for' do professor
# --------------------------------------------------------------------------
n  <- 100000
xt <- 0
for (i in 2:n) {
  xt[i] <- alfa * xt[i - 1] + rnorm(1, 0, sqrt(sigma2))
}
mean(xt)              # deve ficar perto de 0  (E[X] = 0)
var(xt)               # deve ficar perto de 2,78 (Var[X] = sigma^2/(1-alfa^2))

# Mesma coisa com arima.sim (como o professor faz)
ar1 <- arima.sim(n = n, list(ar = c(alfa)))
mean(ar1)
var(ar1)

# --------------------------------------------------------------------------
# 4) E se |alfa| >= 1? A soma explode e a variancia nao existe
# --------------------------------------------------------------------------
sum(1.0^(2 * (0:500)))     # alfa = 1: soma = 501 e continua crescendo
sum(1.1^(2 * (0:100)))     # alfa = 1,1: explode

# FRASE TECNICA: todo AR(1) estacionario (|alfa| < 1) e um MA infinito com
# pesos alfa^j; a variancia e a serie geometrica sigma^2/(1 - alfa^2), que
# so converge porque |alfa| < 1.
