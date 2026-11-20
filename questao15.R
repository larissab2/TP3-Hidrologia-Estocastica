# =====================================================================
# 3o TRABALHO PRATICO - HIDROLOGIA ESTOCASTICA (EHR806)
# QUESTAO 15 - Repetir os itens a) e b) da Questao 14 com a serie
#              transformada  Y = ln(QAH)   (logaritmo natural da vazao anual)
#              (o enunciado indica: ver Salas (1980), pagina 70)
# Salas (1980), p. 70, secao 3.2.1 (Normalizacao de series anuais): varias series
#   anuais de vazao se ajustam bem a uma lognormal de 2 parametros; entao
#   y = log(x) tem distribuicao normal (media mu_y e desvio sigma_y). (So li a p. 70.)
# Baseado em: Aulas-11/03-Selecao de Modelos via Criterio de Informacao.R
#             Aulas-11/04-autoarima.R e 02-ARMA.R ; semente: Aulas/09/ProcEstocasticos.r
# Por que o logaritmo? A vazao e assimetrica (muitos anos "normais" e poucos
# muito cheios) e o AR(1) gaussiano e simetrico. O ln "simetriza" a serie.
# COMO USAR: rode o script (Source). Vao abrir 2 janelinhas: clique no arquivo ANUAL
#   (Qanual.csv) e depois no MENSAL (Qmes.csv). Nao precisa mexer em pasta nenhuma.
# =====================================================================
library(tseries)
library(forecast)
graphics.off()
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
q   <- anual$Qanual                  # vazao anual original
ano <- anual$ANO
n   <- length(q)
y   <- log(q)                        # Y = ln(QAH)  (log() no R ja e o logaritmo natural)

# Funcao pequena para o coeficiente de assimetria (Salas)
assim <- function(z) { nn <- length(z); nn/((nn - 1)*(nn - 2)) * sum(((z - mean(z))/sd(z))^3) }

# =====================================================================
# a) ESTATISTICAS DESCRITIVAS E GRAFICOS DE Y
# =====================================================================
my <- mean(y); sy <- sd(y)
desc <- rbind(
  original = c(media = mean(q), desvio = sd(q), assimetria = assim(q), minimo = min(q), maximo = max(q)),
  ln.Y     = c(media = my,      desvio = sy,    assimetria = assim(y), minimo = min(y), maximo = max(y)))
print(round(desc, 3))
summary(y)
cat("Autocorrelacao lag 1: Q =", round(acf(q, plot = FALSE)$acf[2], 3),
    "| Y =", round(acf(y, plot = FALSE)$acf[2], 3), "\n")

# Teste de normalidade: antes e depois do logaritmo (H0 = normal)
shapiro.test(q)
shapiro.test(y)

# Estatisticas mensais do ln da vazao mensal
mm <- as.matrix(mes[, -1])
lmm <- log(mm)
est.mes <- data.frame(media = colMeans(lmm), desvio = apply(lmm, 2, sd),
                      minimo = apply(lmm, 2, min), maximo = apply(lmm, 2, max))
print(round(est.mes, 3))

# Grafico 1: serie anual de Y
plot(ano, y, type = "l", col = "red", xlab = "ano", ylab = "Y = ln(vazao anual)",
     main = "Q15a - Serie anual transformada Y = ln(QAH)")
abline(h = my, col = "blue", lty = 2)

# Grafico 2: serie mensal de ln(vazao)
qm.ts <- ts(as.vector(t(lmm)), start = c(ano[1], 1), frequency = 12)
plot(qm.ts, col = "red", xlab = "ano", ylab = "ln(vazao mensal)",
     main = "Q15a - ln das vazoes medias mensais (1938-2013)")

# Grafico 3: media mensal de ln(vazao) +- 1 desvio
meses <- c("JAN","FEV","MAR","ABR","MAI","JUN","JUL","AGO","SET","OUT","NOV","DEZ")
plot(1:12, est.mes$media, type = "b", pch = 19, col = "blue", xaxt = "n",
     ylim = range(est.mes$media - est.mes$desvio, est.mes$media + est.mes$desvio),
     xlab = "mes", ylab = "ln(vazao)", main = "Q15a - Media mensal de ln(vazao) +- 1 desvio")
axis(1, at = 1:12, labels = meses)
arrows(1:12, est.mes$media - est.mes$desvio, 1:12, est.mes$media + est.mes$desvio,
       angle = 90, code = 3, length = 0.05, col = "gray40")

# Grafico 4: o efeito do logaritmo (histogramas antes e depois)
par(mfrow = c(1, 2))
hist(q, col = "lightblue", main = "Q (original): assimetrica", xlab = "vazao anual")
hist(y, col = "lightgreen", main = "Y = ln(Q): quase simetrica", xlab = "ln(vazao anual)")
par(mfrow = c(1, 1))

# =====================================================================
# b) SELECAO DO MODELO PARA Y
# =====================================================================
kpss.test(y, null = "L")             # H0 = estacionaria
adf.test(y)                          # H0 = raiz unitaria

par(mfrow = c(1, 3))
plot(ano, y, type = "l", col = "red", main = "Y = ln(Q)", xlab = "ano", ylab = "Y")
acf(y, type = "correlation", lag.max = 20, main = "FAC (ACF)")
pacf(y, lag.max = 20, main = "FACP (PACF)")
par(mfrow = c(1, 1))

x <- y - my                          # serie centralizada (como o professor faz com o lynx)

# AR(p), p = 1 a 5 (formulas do professor)
p <- 1:5
crit <- matrix(0, nrow = length(p), ncol = 3)
for (k in p) {
  ar <- arima(x, order = c(k, 0, 0), include.mean = FALSE)
  crit[k, 1] <- k
  crit[k, 2] <- n*log(ar$sigma2) + 2*k
  crit[k, 3] <- n*log(ar$sigma2) + log(n)*k
}
colnames(crit) <- c("p", "AIC", "BIC")
print(round(crit, 2))

# ARMA(p,q), p = 0..3 e q = 0..3
tab <- NULL
for (pp in 0:3) {
  for (qq in 0:3) {
    f <- try(arima(x, order = c(pp, 0, qq), include.mean = FALSE), silent = TRUE)
    if (inherits(f, "try-error")) next
    k <- pp + qq
    K <- k + 1
    tab <- rbind(tab, data.frame(p = pp, q = qq,
            AIC  = n*log(f$sigma2) + 2*k,
            BIC  = n*log(f$sigma2) + log(n)*k,
            AICc = AIC(f) + 2*K*(K + 1)/(n - K - 1)))   # outra escala; vale a ORDEM
  }
}
tab <- tab[order(tab$AIC), ]
print(round(tab, 2))
cat("Menor AIC : (", tab$p[which.min(tab$AIC)], ",", tab$q[which.min(tab$AIC)], ")\n")
cat("Menor BIC : (", tab$p[which.min(tab$BIC)], ",", tab$q[which.min(tab$BIC)], ")\n")
cat("Menor AICc: (", tab$p[which.min(tab$AICc)], ",", tab$q[which.min(tab$AICc)], ")\n")

# Conferencia com o auto.arima (como no 04-autoarima.R)
mod.auto <- auto.arima(y, trace = TRUE)
print(mod.auto)

# Modelo escolhido: AR(1)
mod <- arima(x, order = c(1, 0, 0), include.mean = FALSE)
print(mod)
alfa   <- as.numeric(mod$coef["ar1"])
sigma2 <- mod$sigma2
cat("alfa1 =", round(alfa, 4), " | sigma2 =", round(sigma2, 5),
    " | erro padrao =", round(sqrt(mod$var.coef[1, 1]), 4), "\n")

# Validacao dos residuos
checkresiduals(mod)
kpss.test(mod$residuals, null = "L")
mean(mod$residuals)
sd(mod$residuals)
shapiro.test(mod$residuals)          # agora os residuos SAO normais (p alto)
par(mfrow = c(1, 2))
acf(mod$residuals, lag.max = 20, main = "FAC dos residuos")
pacf(mod$residuals, lag.max = 20, main = "FACP dos residuos")
par(mfrow = c(1, 1))

# =====================================================================
# b) GERACAO DE 100 SERIES DE 76 ANOS
# =====================================================================
# 1) gera Y (media + AR(1) simulado)   2) volta para vazao com exp()
nser <- 100
Yger <- matrix(0, nrow = n, ncol = nser)
for (i in 1:nser) {
  xg <- arima.sim(n = n, list(ar = alfa), sd = sqrt(sigma2))
  Yger[, i] <- my + as.numeric(xg)               # serie gerada de Y
}
Qger <- exp(Yger)                                 # volta para a escala original (m3/s)
dim(Qger)

comp <- data.frame(
  estatistica = c("media", "desvio padrao", "assimetria", "autocorr. lag 1", "minimo", "maximo"),
  observada = c(mean(q), sd(q), assim(q), acf(q, plot = FALSE)$acf[2], min(q), max(q)),
  geradas.media.das.100 = c(mean(colMeans(Qger)), mean(apply(Qger, 2, sd)),
        mean(apply(Qger, 2, assim)), mean(apply(Qger, 2, function(z) acf(z, plot = FALSE)$acf[2])),
        mean(apply(Qger, 2, min)), mean(apply(Qger, 2, max))))
cat("--- escala ORIGINAL (m3/s) ---\n"); print(comp, digits = 4)
compY <- data.frame(
  estatistica = c("media", "desvio padrao", "autocorr. lag 1"),
  observada = c(my, sy, acf(y, plot = FALSE)$acf[2]),
  geradas.media.das.100 = c(mean(colMeans(Yger)), mean(apply(Yger, 2, sd)),
        mean(apply(Yger, 2, function(z) acf(z, plot = FALSE)$acf[2]))))
cat("--- escala LOGARITMICA (Y) ---\n"); print(compY, digits = 4)
cat("Series geradas com algum valor negativo:", sum(apply(Qger, 2, min) < 0), "(com exp() e impossivel)\n")

matplot(ano, Qger, type = "l", lty = 1, col = "gray80", xlab = "ano", ylab = "vazao (m3/s)",
        main = "Q15b - 100 series geradas (cinza) e observada (vermelho)")
lines(ano, q, col = "red", lwd = 2)

par(mfrow = c(1, 2))
hist(colMeans(Qger), col = "lightblue", main = "Medias das 100 series", xlab = "media")
abline(v = mean(q), col = "red", lwd = 2)
hist(apply(Qger, 2, sd), col = "lightblue", main = "Desvios das 100 series", xlab = "desvio padrao")
abline(v = sd(q), col = "red", lwd = 2)
par(mfrow = c(1, 1))

colnames(Qger) <- paste0("S", 1:nser)
write.table(data.frame(ANO = ano, Qger), "series_geradas_Q15.csv",
            sep = ";", dec = ".", row.names = FALSE)

# FRASE TECNICA: a transformacao Y = ln(Q) reduz a assimetria (0,70 -> 0,04) e torna a
# serie compativel com a normal (Shapiro p = 0,94 contra 0,025 da serie original). Entre
# AR(p) e ARMA(p,q), AIC, BIC e auto.arima escolhem o AR(1) para Y (alfa1 = 0,34), com
# residuos ruido branco (Ljung-Box p = 0,58) e normais (Shapiro p = 0,52). As 100 series
# geradas em Y e levadas de volta por exp() reproduzem media, desvio, assimetria e
# limites da serie observada e nao produzem vazoes negativas.
