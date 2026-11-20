# =====================================================================
# 3o TRABALHO PRATICO - HIDROLOGIA ESTOCASTICA (EHR806)
# QUESTAO 14 - Serie de vazoes (76 anos, 1938-2013)
#   a) Estatisticas descritivas e graficos da serie anual e das medias mensais
#   b) Selecionar (de forma justificada) entre AR(p) e ARMA(p,q), ajustar na serie
#      anual e gerar 100 series sinteticas de 76 anos
# =====================================================================
library(tseries)
library(forecast)
graphics.off()                       # fecha graficos antigos
set.seed(13679821)                   # semente do professor

# ---------------------------------------------------------------------
# ESCOLHER OS ARQUIVOS (uma janelinha abre para voce clicar no arquivo)
# ---------------------------------------------------------------------
cat("\n>>> Escolha o arquivo ANUAL (o que tem a coluna Qanual) na janela que vai abrir...\n")
arq_anual <- file.choose()
cat(">>> Agora escolha o arquivo MENSAL (o que tem JAN, FEV, MAR ...) na janela que vai abrir...\n")
arq_mes   <- file.choose()
setwd(dirname(arq_anual))            # as series geradas serao salvas na pasta do arquivo anual

# ---------------------------------------------------------------------
# LEITURA DOS DADOS  (separador ";" e decimal com PONTO)
# ---------------------------------------------------------------------
anual <- read.table(arq_anual, header = TRUE, sep = ";", dec = ".")
mes   <- read.table(arq_mes,   header = TRUE, sep = ";", dec = ".")
if (!("Qanual" %in% names(anual)))   # confere se nao trocou os arquivos
  stop("O arquivo ANUAL escolhido nao tem a coluna Qanual. Voce escolheu o arquivo certo? Rode de novo.")
if (!("JAN" %in% names(mes)))
  stop("O arquivo MENSAL escolhido nao tem a coluna JAN. Voce escolheu o arquivo certo? Rode de novo.")
q   <- anual$Qanual                  # vazoes medias anuais (m3/s)
ano <- anual$ANO
n   <- length(q)                     # 76 anos
str(anual)

# =====================================================================
# a) ESTATISTICAS DESCRITIVAS E GRAFICOS
# =====================================================================
m  <- mean(q)
s  <- sd(q)
g  <- n/((n - 1)*(n - 2)) * sum(((q - m)/s)^3)            # coef. de assimetria (Salas)
desc <- c(n = n, minimo = min(q), maximo = max(q), media = m, mediana = median(q),
          desvio.padrao = s, coef.variacao = s/m, assimetria = g,
          autocorr.lag1 = acf(q, plot = FALSE)$acf[2])
print(round(desc, 3))
summary(q)

# Estatisticas por mes (media, desvio padrao, minimo, maximo) -- 76 valores por mes
mm <- as.matrix(mes[, -1])                                  # tira a coluna ANO
est.mes <- data.frame(media = colMeans(mm), desvio = apply(mm, 2, sd),
                      minimo = apply(mm, 2, min), maximo = apply(mm, 2, max))
print(round(est.mes, 2))

# Grafico 1: serie anual
plot(ano, q, type = "l", col = "red", xlab = "ano", ylab = "vazao media anual (m3/s)",
     main = "Q14a - Vazoes medias anuais (1938-2013)")
abline(h = m, col = "blue", lty = 2)            # linha da media

# Grafico 2: serie mensal completa (76 anos x 12 meses = 912 valores)
qm <- as.vector(t(mm))                          # coloca mes apos mes, ano apos ano
qm.ts <- ts(qm, start = c(ano[1], 1), frequency = 12)
plot(qm.ts, col = "red", xlab = "ano", ylab = "vazao media mensal (m3/s)",
     main = "Q14a - Vazoes medias mensais (1938-2013)")

# Grafico 3: media mensal (o "ciclo do ano") com +- 1 desvio padrao
meses <- c("JAN","FEV","MAR","ABR","MAI","JUN","JUL","AGO","SET","OUT","NOV","DEZ")
plot(1:12, est.mes$media, type = "b", pch = 19, col = "blue", xaxt = "n",
     ylim = range(est.mes$media - est.mes$desvio, est.mes$media + est.mes$desvio),
     xlab = "mes", ylab = "vazao (m3/s)", main = "Q14a - Media mensal +- 1 desvio padrao")
axis(1, at = 1:12, labels = meses)
arrows(1:12, est.mes$media - est.mes$desvio, 1:12, est.mes$media + est.mes$desvio,
       angle = 90, code = 3, length = 0.05, col = "gray40")

# =====================================================================
# b) SELECAO DO MODELO (serie anual)
# =====================================================================
# --- 1) A serie e estacionaria? (KPSS: H0 = estacionaria ; ADF: H0 = raiz unitaria)
kpss.test(q, null = "L")
adf.test(q)
shapiro.test(q)                      # a serie e normal? (H0 = normal)

# --- 2) FAC e FACP da serie (para ter uma ideia do modelo)
par(mfrow = c(1, 3))
plot(ano, q, type = "l", col = "red", main = "Serie anual", xlab = "ano", ylab = "m3/s")
acf(q, type = "correlation", lag.max = 20, main = "FAC (ACF)")
pacf(q, lag.max = 20, main = "FACP (PACF)")
par(mfrow = c(1, 1))

# --- 3) Centralizar a serie (como o professor faz com o lynx): x = q - media
x <- q - m

# --- 4) AR(p), p = 1 a 5: AIC e BIC com as formulas do professor
p <- 1:5
crit <- matrix(0, nrow = length(p), ncol = 3)
for (k in p) {
  ar <- arima(x, order = c(k, 0, 0), include.mean = FALSE)
  crit[k, 1] <- k
  crit[k, 2] <- n*log(ar$sigma2) + 2*k              # AIC do professor
  crit[k, 3] <- n*log(ar$sigma2) + log(n)*k         # BIC do professor
}
colnames(crit) <- c("p", "AIC", "BIC")
print(round(crit, 2))

# --- 5) ARMA(p,q): testar todas as combinacoes p = 0..3 e q = 0..3
tab <- NULL
for (pp in 0:3) {
  for (qq in 0:3) {
    f <- try(arima(x, order = c(pp, 0, qq), include.mean = FALSE), silent = TRUE)
    if (inherits(f, "try-error")) next                # se o R nao conseguiu ajustar, pula
    k <- pp + qq
    K <- k + 1                                        # parametros contando a variancia
    tab <- rbind(tab, data.frame(p = pp, q = qq,
            AIC  = n*log(f$sigma2) + 2*k,
            BIC  = n*log(f$sigma2) + log(n)*k,
            AICc = AIC(f) + 2*K*(K + 1)/(n - K - 1)))  # AICc (slides da Aula 11); vem do AIC() do R,
                                                       # por isso esta em outra escala; o que vale e a ORDEM
  }
}
tab <- tab[order(tab$AIC), ]                          # do melhor para o pior
print(round(tab, 2))
cat("Menor AIC : AR/ARMA (", tab$p[which.min(tab$AIC)], ",", tab$q[which.min(tab$AIC)], ")\n")
cat("Menor BIC : AR/ARMA (", tab$p[which.min(tab$BIC)], ",", tab$q[which.min(tab$BIC)], ")\n")
cat("Menor AICc: AR/ARMA (", tab$p[which.min(tab$AICc)], ",", tab$q[which.min(tab$AICc)], ")\n")

# --- 6) Conferencia: o auto.arima escolhe a mesma coisa? (como no 04-autoarima.R)
mod.auto <- auto.arima(q, trace = TRUE)
print(mod.auto)
mod.auto2 <- auto.arima(q, trace = FALSE, stepwise = FALSE, approximation = FALSE)
print(mod.auto2)

# --- 7) MODELO ESCOLHIDO: pegamos o de menor criterio (AR(1))
mod <- arima(x, order = c(1, 0, 0), include.mean = FALSE)
print(mod)
alfa   <- as.numeric(mod$coef["ar1"])
sigma2 <- mod$sigma2
cat("alfa1 =", round(alfa, 4), " | sigma2 =", round(sigma2, 2),
    " | erro padrao de alfa1 =", round(sqrt(mod$var.coef[1, 1]), 4), "\n")
cat("Raiz de 1 - alfa*B fora do circulo unitario? |1/alfa| =", round(abs(1/alfa), 2), "\n")

# --- 8) VALIDACAO: os residuos parecem ruido branco? (como no 04-autoarima.R)
checkresiduals(mod)                  # grafico + teste de Ljung-Box
kpss.test(mod$residuals, null = "L")
mean(mod$residuals)
sd(mod$residuals)
shapiro.test(mod$residuals)
par(mfrow = c(1, 2))
acf(mod$residuals, lag.max = 20, main = "FAC dos residuos")
pacf(mod$residuals, lag.max = 20, main = "FACP dos residuos")
par(mfrow = c(1, 1))

# =====================================================================
# b) GERACAO DE 100 SERIES SINTETICAS DE 76 ANOS
# =====================================================================
# Cada serie: Qger = media + AR(1) simulado   (volta a media que tiramos no passo 3)
nser <- 100
Qger <- matrix(0, nrow = n, ncol = nser)
for (i in 1:nser) {
  xg <- arima.sim(n = n, list(ar = alfa), sd = sqrt(sigma2))   # como no 02-ARMA.R
  Qger[, i] <- m + as.numeric(xg)
}
dim(Qger)                                            # 76 linhas (anos) x 100 colunas (series)

# Comparar as series geradas com a observada
comp <- data.frame(
  estatistica = c("media", "desvio padrao", "assimetria", "autocorr. lag 1", "minimo", "maximo"),
  observada = c(m, s, g, acf(q, plot = FALSE)$acf[2], min(q), max(q)),
  geradas.media.das.100 = c(mean(colMeans(Qger)), mean(apply(Qger, 2, sd)),
        mean(apply(Qger, 2, function(z) n/((n-1)*(n-2))*sum(((z-mean(z))/sd(z))^3))),
        mean(apply(Qger, 2, function(z) acf(z, plot = FALSE)$acf[2])),
        mean(apply(Qger, 2, min)), mean(apply(Qger, 2, max))))
print(comp, digits = 4)
cat("Series geradas com algum valor negativo:", sum(apply(Qger, 2, min) < 0), "de", nser, "\n")

# Grafico: as 100 series (cinza) e a observada (vermelho)
matplot(ano, Qger, type = "l", lty = 1, col = "gray80", xlab = "ano", ylab = "vazao (m3/s)",
        main = "Q14b - 100 series geradas (cinza) e a serie observada (vermelho)")
lines(ano, q, col = "red", lwd = 2)

# Grafico: a media e o desvio de cada serie gerada contra os da observada
par(mfrow = c(1, 2))
hist(colMeans(Qger), col = "lightblue", main = "Medias das 100 series", xlab = "media")
abline(v = m, col = "red", lwd = 2)
hist(apply(Qger, 2, sd), col = "lightblue", main = "Desvios das 100 series", xlab = "desvio padrao")
abline(v = s, col = "red", lwd = 2)
par(mfrow = c(1, 1))

# Guardar as 100 series em arquivo (colunas S1 ... S100, uma linha por ano)
colnames(Qger) <- paste0("S", 1:nser)
write.table(data.frame(ANO = ano, Qger), "series_geradas_Q14.csv",
            sep = ";", dec = ".", row.names = FALSE)
