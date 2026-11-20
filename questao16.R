# =====================================================================
# 3o TRABALHO PRATICO - HIDROLOGIA ESTOCASTICA (EHR806)
# QUESTAO 16 - Serie MENSAL (76 anos = 912 meses)
#   - selecionar (de forma justificada) um modelo SARIMA(p,d,q)(P,D,Q)12 e ajusta-lo
#   - gerar 100 series de 76 anos
#   - comparar MEDIAS MENSAIS e DESVIOS PADRAO MENSAIS gerados x observados
# Baseado em: Aulas-11/03-ARMA-ARIMA.pdf (SARIMA: S = 12 para dados mensais;
#               "opcao de remover o elemento sazonal e usar o ARIMA padrao")
#             Aulas-11/04-autoarima.R (auto.arima, checkresiduals, kpss, shapiro)
#             Aulas-11/03-Selecao...R (AIC, BIC) ; semente: Aulas/09/ProcEstocasticos.r
# ESTRUTURA:  PARTE 1 = SARIMA (o que o enunciado pede)
#             PARTE 2 = alternativa do proprio slide: tirar a sazonalidade (padronizar
#                       cada mes) e usar ARMA comum  (EXTRA, para comparar)
# COMO USAR: rode o script (Source). Vai abrir uma janelinha: clique no arquivo MENSAL
#   (Qmes.csv). Nao precisa mexer em pasta nenhuma.
# ATENCAO: o script leva alguns minutos (a grade de modelos SARIMA e demorada).
# =====================================================================
library(tseries)
library(forecast)
graphics.off()
set.seed(13679821)                   # semente do professor

# ---------------------------------------------------------------------
# ESCOLHER O ARQUIVO (uma janelinha abre para voce clicar no arquivo)
# ---------------------------------------------------------------------
cat("\n>>> Escolha o arquivo MENSAL (o que tem JAN, FEV, MAR ...) na janela que vai abrir...\n")
arq_mes <- file.choose()
setwd(dirname(arq_mes))              # as series geradas serao salvas na pasta do arquivo

mes <- read.table(arq_mes, header = TRUE, sep = ";", dec = ".")
if (!("JAN" %in% names(mes)))        # confere se nao escolheu o arquivo errado
  stop("O arquivo MENSAL escolhido nao tem a coluna JAN. Voce escolheu o arquivo certo? Rode de novo.")
mm  <- as.matrix(mes[, -1])          # 76 linhas (anos) x 12 colunas (meses)
ano <- mes$ANO
na  <- nrow(mm)                      # 76 anos
nmes <- na*12                        # 912 meses
meses <- c("JAN","FEV","MAR","ABR","MAI","JUN","JUL","AGO","SET","OUT","NOV","DEZ")

# A serie mensal como UMA fila de 912 numeros: janeiro/1938, fev/1938, ..., dez/2013
z <- ts(as.vector(t(mm)), start = c(ano[1], 1), frequency = 12)

# Estatisticas mensais OBSERVADAS (e o que vamos tentar reproduzir)
obs.media  <- colMeans(mm)
obs.desvio <- apply(mm, 2, sd)
print(round(rbind(media = obs.media, desvio = obs.desvio), 2))

# Funcao que compara as series geradas com a observada (mes a mes)
# 'ger' e uma lista com 100 matrizes 76 x 12
compara <- function(ger, titulo) {
  M <- t(sapply(ger, colMeans))                  # 100 x 12 : media de cada mes em cada serie
  S <- t(sapply(ger, function(g) apply(g, 2, sd)))
  tab <- data.frame(mes = meses, media.obs = obs.media, media.ger = colMeans(M),
                    erro.media.pct = 100*(colMeans(M) - obs.media)/obs.media,
                    desvio.obs = obs.desvio, desvio.ger = colMeans(S),
                    erro.desvio.pct = 100*(colMeans(S) - obs.desvio)/obs.desvio)
  cat("\n=====", titulo, "=====\n"); print(tab, digits = 4, row.names = FALSE)
  cat("Erro relativo medio (valor absoluto):  medias =",
      round(mean(abs(tab$erro.media.pct)), 1), "%   desvios =",
      round(mean(abs(tab$erro.desvio.pct)), 1), "%\n")
  cat("Series com algum valor negativo:", sum(sapply(ger, function(g) any(g < 0))), "de", length(ger), "\n")
  par(mfrow = c(1, 2))
  plot(1:12, obs.media, type = "b", pch = 19, col = "red", lwd = 2, xaxt = "n",
       ylim = range(M, obs.media), xlab = "mes", ylab = "vazao media (m3/s)",
       main = paste0(titulo, "\nmedias mensais"), cex.main = 1)
  axis(1, at = 1:12, labels = meses)
  for (i in 1:nrow(M)) lines(1:12, M[i, ], col = "gray80")
  lines(1:12, obs.media, type = "b", pch = 19, col = "red", lwd = 2)
  legend("top", legend = c("observada", "100 series geradas"), col = c("red", "gray60"),
         lty = 1, lwd = c(2, 1), pch = c(19, NA), bty = "n", cex = 0.8)
  plot(1:12, obs.desvio, type = "b", pch = 19, col = "red", lwd = 2, xaxt = "n",
       ylim = range(S, obs.desvio), xlab = "mes", ylab = "desvio padrao (m3/s)",
       main = paste0(titulo, "\ndesvios mensais"), cex.main = 1)
  axis(1, at = 1:12, labels = meses)
  for (i in 1:nrow(S)) lines(1:12, S[i, ], col = "gray80")
  lines(1:12, obs.desvio, type = "b", pch = 19, col = "red", lwd = 2)
  par(mfrow = c(1, 1))
  invisible(tab)
}

# =====================================================================
# PARTE 1 - SARIMA(p,d,q)(P,D,Q)12
# =====================================================================
# --- 1) FAC e FACP da serie mensal: aparece a sazonalidade? (picos nos lags 12, 24, 36)
par(mfrow = c(1, 2))
acf(z, lag.max = 48, main = "FAC da serie mensal")
pacf(z, lag.max = 48, main = "FACP da serie mensal")
par(mfrow = c(1, 1))

# --- 2) Quantas diferencas? (d = diferenca simples ; D = diferenca sazonal)
kpss.test(as.numeric(z), null = "L")       # H0 = estacionaria
adf.test(as.numeric(z))                    # H0 = raiz unitaria
cat("ndiffs (d) =", ndiffs(z), "  nsdiffs (D) =", nsdiffs(z), "\n")
# d = 0 (nao precisa diferenca simples) e D = 1 (precisa uma diferenca sazonal)

# --- 3) auto.arima (como no 04-autoarima.R)
mod.auto <- auto.arima(z, trace = FALSE)
print(mod.auto)
# CUIDADO: a busca "stepwise" (passo a passo) pode parar num modelo que nao e o
# melhor. Por isso, abaixo, testamos TODAS as combinacoes de uma grade.
# (opcional, demora varios minutos:)
# mod.auto2 <- auto.arima(z, trace = TRUE, stepwise = FALSE, approximation = FALSE)

# --- 4) GRADE de modelos: p,q = 0..2 ; P,Q = 0..1 ; d = 0 ; D = 0 ou 1
tab <- NULL
for (D in 0:1) for (p in 0:2) for (q in 0:2) for (P in 0:1) for (Q in 0:1) {
  f <- try(suppressWarnings(Arima(z, order = c(p, 0, q),
           seasonal = list(order = c(P, D, Q), period = 12), include.mean = (D == 0))),
           silent = TRUE)
  if (inherits(f, "try-error")) next                  # modelo que o R nao consegue ajustar
  tab <- rbind(tab, data.frame(p = p, d = 0, q = q, P = P, D = D, Q = Q,
                               AIC = f$aic, AICc = f$aicc, BIC = f$bic))
}
tab <- tab[order(tab$AICc), ]
print(head(tab, 10), row.names = FALSE)
b <- tab[1, ]
cat("MODELO ESCOLHIDO (menor AICc): SARIMA(", b$p, ",", b$d, ",", b$q, ")(", b$P, ",", b$D, ",", b$Q, ")12\n")
cat("Menor BIC escolhe: SARIMA(", tab$p[which.min(tab$BIC)], ",0,", tab$q[which.min(tab$BIC)], ")(",
    tab$P[which.min(tab$BIC)], ",", tab$D[which.min(tab$BIC)], ",", tab$Q[which.min(tab$BIC)], ")12\n")

# --- 5) AJUSTE do modelo escolhido
mod <- Arima(z, order = c(b$p, 0, b$q), seasonal = list(order = c(b$P, b$D, b$Q), period = 12),
             include.mean = (b$D == 0))
print(mod)

# --- 6) VALIDACAO: os residuos parecem ruido branco? (como no 04-autoarima.R)
checkresiduals(mod)
kpss.test(as.numeric(mod$residuals), null = "L")
mean(mod$residuals)
sd(mod$residuals)
shapiro.test(as.numeric(mod$residuals))

# --- 7) GERAR 100 series de 76 anos (912 meses) com o modelo ajustado
# O arima.sim() do R nao faz diferenca SAZONAL, entao montamos a serie em 2 passos:
#   (a) gera a parte ARMA "w" (com os polinomios sazonais multiplicados)
#   (b) "desfaz" a diferenca sazonal:  x(t) = x(t-12) + w(t),
#       partindo dos 12 primeiros meses OBSERVADOS (1o ano) como valores iniciais.
# (Obs.: testei o simulate() do pacote forecast; ele parte do zero e perde o ciclo
#  sazonal das medias; por isso preferi montar a equacao do modelo na mao.)
cf <- coef(mod)
phi   <- as.numeric(cf[grep("^ar",  names(cf))])      # AR comum
theta <- as.numeric(cf[grep("^ma",  names(cf))])      # MA comum
Phi   <- as.numeric(cf[grep("^sar", names(cf))])      # AR sazonal
Theta <- as.numeric(cf[grep("^sma", names(cf))])      # MA sazonal
multiplica <- function(a, b) convolve(a, rev(b), type = "open")   # multiplica dois polinomios
sazonal <- function(coef, sinal) {                    # monta 1 +- c1 B^12 +- c2 B^24 ...
  v <- c(1, rep(0, 12*length(coef)))
  for (k in seq_along(coef)) v[12*k + 1] <- sinal*coef[k]
  v
}
pol.ar <- multiplica(c(1, -phi),  sazonal(Phi, -1))   # (1 - phi B ...)(1 - Phi B^12 ...)
pol.ma <- multiplica(c(1, theta), sazonal(Theta, +1)) # (1 + theta B ...)(1 + Theta B^12 ...)
modelo.sim <- list()
if (length(pol.ar) > 1) modelo.sim$ar <- -pol.ar[-1]
if (length(pol.ma) > 1) modelo.sim$ma <-  pol.ma[-1]

nser <- 100
ger1 <- vector("list", nser)
for (i in 1:nser) {
  w <- as.numeric(arima.sim(modelo.sim, n = nmes, sd = sqrt(mod$sigma2)))
  x <- w
  if (b$D == 1) {                                     # desfaz a diferenca sazonal
    x <- numeric(nmes)
    x[1:12] <- mm[1, ]                                # 1o ano observado
    for (t in 13:nmes) x[t] <- x[t - 12] + w[t]
  } else {
    x <- x + as.numeric(cf["intercept"])              # sem diferenca: soma a media
  }
  ger1[[i]] <- matrix(x, ncol = 12, byrow = TRUE)     # 76 anos x 12 meses
}

# --- 8) COMPARAR medias e desvios mensais
tab1 <- compara(ger1, "Q16 SARIMA")

# =====================================================================
# PARTE 2 (EXTRA) - tirar a sazonalidade e usar ARMA comum (opcao do slide)
#   2.1  Y = ln(vazao)       -> a vazao mensal e muito assimetrica (como na Q15)
#   2.2  padroniza CADA MES:  w = (Y - media do mes) / desvio do mes
#        (assim cada mes passa a ter media 0 e desvio 1: a sazonalidade some)
#   2.3  ajusta um ARMA(p,q) em w  (escolha por AICc)
#   2.4  gera w, "volta" (multiplica pelo desvio, soma a media do mes) e aplica exp()
# =====================================================================
Y   <- log(mm)
muY <- colMeans(Y)
sgY <- apply(Y, 2, sd)
w   <- as.vector(t(sweep(sweep(Y, 2, muY), 2, sgY, "/")))    # serie padronizada (912 valores)
par(mfrow = c(1, 2)); acf(w, lag.max = 48, main = "FAC de w (sem sazonalidade)")
pacf(w, lag.max = 48, main = "FACP de w"); par(mfrow = c(1, 1))
kpss.test(w, null = "L")

tab2 <- NULL
for (p in 0:3) for (q in 0:3) {
  f <- try(Arima(w, order = c(p, 0, q), include.mean = FALSE), silent = TRUE)
  if (inherits(f, "try-error")) next
  tab2 <- rbind(tab2, data.frame(p = p, q = q, AIC = f$aic, AICc = f$aicc, BIC = f$bic))
}
tab2 <- tab2[order(tab2$AICc), ]
print(head(tab2, 8), row.names = FALSE)
b2 <- tab2[1, ]
cat("ARMA escolhido (menor AICc): (", b2$p, ",", b2$q, ")\n")
mod2 <- Arima(w, order = c(b2$p, 0, b2$q), include.mean = FALSE)
print(mod2)
checkresiduals(mod2)
shapiro.test(as.numeric(mod2$residuals))

ger2 <- vector("list", nser)
for (i in 1:nser) {
  wg <- as.numeric(simulate(mod2, nsim = nmes, future = FALSE))
  g  <- matrix(wg, ncol = 12, byrow = TRUE)
  g  <- sweep(sweep(g, 2, sgY, "*"), 2, muY, "+")             # volta a escala de Y
  ger2[[i]] <- exp(g)                                          # volta a vazao (m3/s)
}
tab2c <- compara(ger2, "Q16 ARMA padronizado")

# Guardar as series geradas: uma coluna por serie (912 meses), das duas partes
saida1 <- sapply(ger1, function(g) as.vector(t(g)))
saida2 <- sapply(ger2, function(g) as.vector(t(g)))
colnames(saida1) <- colnames(saida2) <- paste0("S", 1:nser)
write.table(saida1, "series_geradas_Q16_SARIMA.csv", sep = ";", dec = ".", row.names = FALSE)
write.table(saida2, "series_geradas_Q16_ARMA_padronizado.csv", sep = ";", dec = ".", row.names = FALSE)

# FRASE TECNICA: a serie mensal (912 valores) tem forte sazonalidade (picos da FAC nos
# lags 12, 24, ...) e pede d = 0 e D = 1 (diferenca sazonal). Na grade, o menor AICc escolhe
# o SARIMA(2,0,2)(0,1,1)12 (o BIC escolhe o (1,0,1)(0,1,1)12), com residuos descorrelacionados
# (Ljung-Box p = 0,92) mas nao normais. Porem, as 100 series geradas reproduzem mal as
# estatisticas mensais: erro medio de 35% nas medias e de 111% nos desvios (todos os meses
# ficam com desvio ~40 m3/s, pois o modelo tem uma unica variancia de inovacao, e as 100
# series tem vazoes negativas). A alternativa do proprio slide (padronizar cada mes e usar um
# ARMA(3,1) no ln da vazao) reproduz medias (erro 0,5%) e desvios (3,5%) sem vazoes negativas.
