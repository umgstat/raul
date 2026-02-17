


# ==========================================================
# analysis.R — RAUL reproducible analysis
# R version 4.4.1 (2024-06-14) -- "Race for Your Life"
# Copyright (C) 2024 The R Foundation for Statistical Computing
# Platform: aarch64-apple-darwin20
# R Studio
# Version 2026.01.0+392 (2026.01.0+392)
# ==========================================================

# ------------------------------
# SECTION 1 — Setup
# ------------------------------
set.seed(123)  # reproducibility

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(pROC)
  library(psych)
})

# Output dirs
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("tables"))  dir.create("tables")

# ------------------------------
# SECTION 2 — Data import
# ------------------------------
url_csv <- "https://raw.githubusercontent.com/umgstat/raul/refs/heads/main/dataset.csv"
rauldataset <- read.csv(url_csv, header = TRUE, stringsAsFactors = TRUE)

# Factor ordering
rauldataset$outcome <- factor(rauldataset$outcome, levels = c("benign", "malignant"))
rauldataset$class1  <- factor(rauldataset$class1,  levels = c("A", "C"))
rauldataset$class2  <- factor(rauldataset$class2,  levels = c("A", "B", "C"))
rauldataset$class3  <- factor(rauldataset$class3,  levels = c("A", "B", "C"))
rauldataset$raul    <- factor(rauldataset$raul,    levels = c("A", "B", "C"))

# Quick integrity check
req_cols <- c("age","outcome","LDHtot","LDH1","LDH2","LDH3","LDH4","LDH5",
              "UMG","classifier1","class1","class2","p","class3","raul")
stopifnot(all(req_cols %in% names(rauldataset)))

# ------------------------------
# SECTION 3 — UMG computation and validation
# ------------------------------
# Check UMG ≈ LDH3 + 24/LDH1
umg_check <- all.equal(rauldataset$UMG, rauldataset$LDH3 + 24/rauldataset$LDH1, tolerance = 1e-2)
print(umg_check)

# Scatter: LDH3 vs 24/LDH1
df_umg <- rauldataset %>%
  mutate(LDH1inv = ifelse(!is.na(LDH1) & LDH1 > 0, 24/LDH1, NA_real_)) %>%
  filter(!is.na(LDH3), !is.na(LDH1inv))

p_ldh3_inv <- ggplot(df_umg, aes(x = LDH3, y = LDH1inv, color = outcome, shape = outcome)) +
  geom_point(alpha = 0.9, size = 2) +
  scale_color_viridis(discrete = TRUE, option = "D", begin = 0.75, end = 0.15) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "LDH3", y = "24 / LDH1", color = "Outcome", shape = "Outcome",
       title = "LDH3 vs 24/LDH1") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())
p_ldh3_inv
# ggsave("figures/ldh3_vs_24overldh1.png", p_ldh3_inv, width = 6.5, height = 4.2, dpi = 300)

# Boxplot: UMG by outcome
p_box <- ggplot(rauldataset, aes(x = outcome, y = UMG, fill = outcome)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.alpha = 0.6) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.75, end = 0.15) +
  labs(x = NULL, y = "UMG index", title = "UMG index by outcome") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")
p_box
# ggsave("figures/box_umg_outcome.png", p_box, width = 6.5, height = 4.2, dpi = 300)

# ------------------------------
# SECTION 4 — Helper: diagnostic metrics
# ------------------------------
metrics_from_cm <- function(cm) {
  TN <- cm["benign","neg"]; FP <- cm["benign","pos"]
  FN <- cm["malignant","neg"]; TP <- cm["malignant","pos"]
  N  <- sum(cm)
  
  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  ppv  <- TP / (TP + FP)
  npv  <- TN / (TN + FN)
  acc  <- (TP + TN) / N
  f1   <- ifelse((2*TP + FP + FN) > 0, 2*TP / (2*TP + FP + FN), NA_real_)
  
  sens_ci <- binom.test(TP, TP + FN)$conf.int
  spec_ci <- binom.test(TN, TN + FP)$conf.int
  ppv_ci  <- binom.test(TP, TP + FP)$conf.int
  npv_ci  <- binom.test(TN, TN + FN)$conf.int
  acc_ci  <- binom.test(TP + TN, N)$conf.int
  
  out <- data.frame(
    Measure  = c("Sensitivity","Specificity","PPV","NPV","Accuracy","F1"),
    Estimate = c(sens, spec, ppv, npv, acc, f1),
    CI_low   = c(sens_ci[1], spec_ci[1], ppv_ci[1], npv_ci[1], acc_ci[1], NA),
    CI_high  = c(sens_ci[2], spec_ci[2], ppv_ci[2], npv_ci[2], acc_ci[2], NA)
  )
  attr(out, "confusion") <- c(TN=TN, FP=FP, FN=FN, TP=TP, N=N)
  return(out)
}


# ------------------------------
# SECTION 5 — Classifier 1 (UMG–LDHtot)
# ------------------------------
# Scatter + boundary (x log10)
p_c1 <- ggplot(rauldataset, aes(x = LDHtot, y = UMG, color = outcome, shape = outcome)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_x_log10() +
  scale_color_viridis(discrete = TRUE, option = "D", begin = 0.75, end = 0.15) +
  scale_shape_manual(values = c(16, 17)) +
  geom_function(fun = function(x) 40 - 0.05*x, color = "black", linewidth = 0.8) +
  labs(x = "total LDH (log10)", y = "UMG index",
       title = "Classifier 1: UMG vs total LDH") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())
p_c1 
# ggsave("figures/classifier1_scatter.png", p_c1, width = 6.5, height = 4.2, dpi = 300)

# Confusion and metrics (A = neg, C = pos)
cm_c1 <- table(rauldataset$outcome, factor(ifelse(rauldataset$class1=="A","neg","pos"),
                                           levels=c("neg","pos")))
print(cm_c1)
perf_c1 <- metrics_from_cm(cm_c1); print(perf_c1)
# write.csv(perf_c1, "tables/classifier1_binary_metrics.csv", row.names = FALSE)

# ROC/AUC on continuous classifier1 (higher = higher risk)
y_bin <- ifelse(rauldataset$outcome=="malignant", 1, 0)
roc_c1 <- roc(response = y_bin, predictor = rauldataset$classifier1, direction = ">")
auc_c1 <- auc(roc_c1); ci_c1 <- ci.auc(roc_c1, method = "delong")
print(auc_c1); print(ci_c1)

png("figures/roc_classifier1.png", width=1600, height=1100, res=200)

plot(roc_c1, col = "#3B82F6", lwd = 3, main = "ROC — Classifier 1 (continuous score)")
abline(a=0, b=1, lty=2, col="gray50")
legend("bottomright",
       legend = paste0("AUC = ", round(as.numeric(auc_c1),3),
                       " (95% CI ", round(as.numeric(ci_c1[1]),3), "–",
                       round(as.numeric(ci_c1[3]),3), ")"),
       bty = "n")

dev.off()

plot(roc_c1, col = "#3B82F6", lwd = 3, main = "ROC — Classifier 1 (continuous score)")
abline(a=0, b=1, lty=2, col="gray50")
legend("bottomright",
       legend = paste0("AUC = ", round(as.numeric(auc_c1),3),
                       " (95% CI ", round(as.numeric(ci_c1[1]),3), "–",
                       round(as.numeric(ci_c1[3]),3), ")"),
       bty = "n")


# ------------------------------
# SECTION 6 — Classifier 2 (UMG–LDH5 triangular rule)
# ------------------------------
df <- rauldataset
xmin <- min(df$UMG, na.rm=TRUE); xmax <- max(df$UMG, na.rm=TRUE)
ymin <- min(df$LDH5, na.rm=TRUE); ymax <- max(df$LDH5, na.rm=TRUE)

rect_green <- data.frame(xmin = xmin, xmax = 29, ymin = ymin, ymax = ymax)
A <- c(29, 1.5); B <- c(29, 13.7); C <- c(32.5, 13.7)
poly_gray <- data.frame(x = c(A[1],B[1],C[1]), y = c(A[2],B[2],C[2]))

p_c2 <- ggplot(df, aes(x = UMG, y = LDH5)) +
  annotate("rect", xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill="red", alpha=0.08) +
  geom_rect(data=rect_green, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            inherit.aes=FALSE, fill="green", alpha=0.08) +
  geom_polygon(data=poly_gray, aes(x=x, y=y), inherit.aes=FALSE, fill="grey60", alpha=0.20) +
  geom_polygon(data=poly_gray, aes(x=x, y=y), inherit.aes=FALSE, fill=NA, color="black", linewidth=0.8) +
  geom_point(aes(color = outcome, shape = outcome), alpha = 0.9, size = 2) +
  scale_color_viridis(discrete = TRUE, option = "D", begin = 0.75, end = 0.15) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "UMG index", y = "LDH5", title = "Classifier 2: UMG–LDH5 risk regions") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())
p_c2
# ggsave("figures/classifier2_triangle.png", p_c2, width = 6.5, height = 4.2, dpi = 300)

# Binary (A=neg; B+C=pos)
cm_c2 <- table(rauldataset$outcome, factor(ifelse(rauldataset$class2=="A","neg","pos"),
                                           levels=c("neg","pos")))
print(cm_c2)
perf_c2 <- metrics_from_cm(cm_c2); print(perf_c2)

# write.csv(perf_c2, "tables/classifier2_binary_metrics.csv", row.names = FALSE)

# ------------------------------
# SECTION 7 — Classifier 3 (simulation-based probability p)
# ------------------------------
# Binary (A=neg; B+C=pos)
cm_c3 <- table(rauldataset$outcome, factor(ifelse(rauldataset$class3=="A","neg","pos"),
                                           levels=c("neg","pos")))
print(cm_c3)
perf_c3 <- metrics_from_cm(cm_c3); print(perf_c3)

# write.csv(perf_c3, "tables/classifier3_binary_metrics.csv", row.names = FALSE)

# ROC/AUC on p (continuous)
y_bin <- ifelse(rauldataset$outcome=="malignant", 0, 1)
roc_p <- roc(response = y_bin, predictor = rauldataset$p, direction = ">")
auc_p <- auc(roc_p); ci_p <- ci.auc(roc_p, method = "delong")
print(auc_p); print(ci_p)

png("figures/roc_classifier3_p.png", width=1600, height=1100, res=200)
plot(roc_p, col = "#7C3AED", lwd = 3, main = "ROC — Classifier 3 (probability p)")
abline(a=0, b=1, lty=2, col="gray50")
legend("bottomright",
       legend = paste0("AUC = ", round(as.numeric(auc_p),3),
                       " (95% CI ", round(as.numeric(ci_p[1]),3), "–",
                       round(as.numeric(ci_p[3]),3), ")"),
       bty = "n")
dev.off()
plot(roc_p, col = "#7C3AED", lwd = 3, main = "ROC — Classifier 3 (probability p)")
abline(a=0, b=1, lty=2, col="gray50")
legend("bottomright",
       legend = paste0("AUC = ", round(as.numeric(auc_p),3),
                       " (95% CI ", round(as.numeric(ci_p[1]),3), "–",
                       round(as.numeric(ci_p[3]),3), ")"),
       bty = "n")

# Append AUC row to binary table
c3_tab <- perf_c3 %>%
  mutate(`95% CI` = ifelse(is.na(CI_low), "",
                           sprintf("[%.3f–%.3f]", CI_low, CI_high))) %>%
  select(Measure, Estimate, `95% CI`)
auc_row3 <- data.frame(Measure="AUC (continuous p)",
                       Estimate=as.numeric(auc_p),
                       `95% CI`=sprintf("[%.3f–%.3f]", as.numeric(ci_p[1]), as.numeric(ci_p[3])))
c3_final <- rbind(c3_tab, auc_row3)
# write.csv(c3_final, "tables/classifier3_perf_with_auc.csv", row.names = FALSE)
print(c3_final)

# ------------------------------
# SECTION 8 — Final RAUL score (A/B/C) and binary A vs B+C
# ------------------------------
# Final table A/B/C
raul_tab <- table(rauldataset$outcome, rauldataset$raul)
print(raul_tab)
# write.csv(as.data.frame(raul_tab), "tables/raul_final_confusion_abc.csv", row.names = FALSE)

# Binary A vs (B+C)
cm_raul <- table(rauldataset$outcome, factor(ifelse(rauldataset$raul=="A","neg","pos"),
                                             levels=c("neg","pos")))
print(cm_raul)
perf_raul <- metrics_from_cm(cm_raul); print(perf_raul)
# write.csv(perf_raul, "tables/raul_binary_metrics.csv", row.names = FALSE)

# Likelihood ratios
TN <- attr(perf_raul, "confusion")["TN"]; FP <- attr(perf_raul, "confusion")["FP"]
FN <- attr(perf_raul, "confusion")["FN"]; TP <- attr(perf_raul, "confusion")["TP"]
sens <- perf_raul$Estimate[perf_raul$Measure=="Sensitivity"]
spec <- perf_raul$Estimate[perf_raul$Measure=="Specificity"]
LR_pos <- sens / (1 - spec)
LR_neg <- (1 - sens) / spec
lr_out <- cbind(LR_pos=LR_pos, LR_neg=LR_neg); print(lr_out)
# write.csv(as.data.frame(lr_out), "tables/raul_likelihood_ratios.csv", row.names = FALSE)

# ------------------------------
# SECTION 9 — Calibration (p) and Agreement (kappa)
# ------------------------------
# Calibration by deciles of p
cal <- rauldataset %>%
  mutate(bin = cut(p, breaks = seq(0,1,by=0.1), include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(mean_p = mean(p, na.rm=TRUE),
            obs_rate = mean(outcome=="malignant", na.rm=TRUE),
            n = n(), .groups="drop")
# cal
# write.csv(cal, "tables/calibration_p_by_decile.csv", row.names = FALSE)

p_cal <- ggplot(cal, aes(x = mean_p, y = obs_rate)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="gray50") +
  geom_point(size=2.6, color="#111827") +
  geom_line(color="#0F766E") +
  labs(x="Mean predicted probability (p)", y="Observed malignancy rate",
       title="Calibration — Classifier 3 (p)") +
  theme_minimal(base_size = 14)
# p_cal
# ggsave("figures/calibration_p.png", p_cal, width = 6.5, height = 4.2, dpi = 300)

# Agreement (Cohen's kappa)
kap12 <- cohen.kappa(cbind(as.numeric(rauldataset$class1),
                           as.numeric(rauldataset$class2)))
kap13 <- cohen.kappa(cbind(as.numeric(rauldataset$class1),
                           as.numeric(rauldataset$class3)))
# print(kap12); print(kap13)

# ------------------------------
# SECTION 10 — Additional figures
# ------------------------------
# Density of UMG by RAUL
p_den <- ggplot(rauldataset, aes(x = UMG, fill = raul)) +
  geom_density(alpha = 0.6) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin=0.1, end=0.9) +
  labs(x="UMG index", y="Density", title="UMG distribution by RAUL class") +
  theme_minimal(base_size = 14)
# p_den
# ggsave("figures/density_umg_by_raul.png", p_den, width = 6.5, height = 4.2, dpi = 300)

# ------------------------------
# SECTION 11 — (Optional) Classifier 3 illustration with simulations
# ------------------------------
# Gaussian benign-like cloud and uniform malignant-like cloud
x_ben <- rauldataset$UMG[rauldataset$outcome=="benign"]
y_ben <- rauldataset$LDH5[rauldataset$outcome=="benign"]
mx <- mean(x_ben); sx <- sd(x_ben)
my <- mean(y_ben); sy <- sd(y_ben)
r  <- cor(x_ben, y_ben)

n <- 1000
z1 <- rnorm(n); z2 <- rnorm(n); z3 <- r*z1 + sqrt(1-r^2)*z2
xB <- mx + sx*z1; yB <- my + sy*z3
xM <- runif(n, 29, 40); yM <- runif(n, 1, 40)

df_sim <- rbind(
  data.frame(x=xB, y=yB, group="green"),
  data.frame(x=xM, y=yM, group="purple")
)

p_sim <- ggplot(df_sim, aes(x=x, y=y, color=group, shape=group)) +
  geom_point(alpha=0.85, size=1.6) +
  scale_color_viridis(discrete=TRUE, option="D", begin=0.75, end=0.15) +
  scale_shape_manual(values=c(16,17)) +
  labs(x="x (simulated UMG)", y="y (simulated LDH5)",
       title="Simulated benign-like (Gaussian) vs malignant-like (Uniform)") +
  theme_minimal(base_size = 14)
p_sim
# ggsave("figures/classifier3_simulated_clouds.png", p_sim, width = 6.5, height = 4.2, dpi = 300)

# Local neighbourhood around a test point
xtest <- 28.5; ytest <- 13.5
df_sim$dist <- sqrt((df_sim$x-xtest)^2 + (df_sim$y-ytest)^2)
k <- ceiling(0.10 * nrow(df_sim))
radius <- sort(df_sim$dist)[k]
df_sim$near <- df_sim$dist < radius

theta <- seq(0, 2*pi, length.out=361)
circle_df <- data.frame(x = xtest + radius*cos(theta),
                        y = ytest + radius*sin(theta))

p_zoom <- ggplot() +
  geom_point(data=df_sim, aes(x=x, y=y), color="grey85", size=1.0, alpha=0.35) +
  geom_point(data = subset(df_sim, near),
             aes(x=x, y=y, color=group, shape=group), alpha=0.9, size=2) +
  geom_path(data=circle_df, aes(x=x, y=y), color="black", linewidth=0.7) +
  annotate("point", x=xtest, y=ytest, shape=4, size=3.5, stroke=1.1, colour="black") +
  coord_cartesian(xlim=c(xtest-radius, xtest+radius),
                  ylim=c(ytest-radius, ytest+radius), expand=FALSE) +
  scale_color_viridis(discrete=TRUE, option="D", begin=0.75, end=0.15) +
  scale_shape_manual(values=c(16,17)) +
  labs(x="x (simulated UMG)", y="y (simulated LDH5)", title="Local neighbourhood (10th percentile radius)") +
  theme_minimal(base_size = 14)
p_zoom
# ggsave("figures/classifier3_local_neighbourhood.png", p_zoom, width = 6.5, height = 4.2, dpi = 300)

# ------------------------------
# SECTION 12 — Session info (for reproducibility logs)
# ------------------------------
# writeLines(capture.output(sessionInfo()), "tables/session_info.txt")









www = "https://raw.githubusercontent.com/umgstat/raul/refs/heads/main/dataset.csv"
rauldataset = read.csv(www, header = TRUE, stringsAsFactors = TRUE)
attach(rauldataset)
head(rauldataset)

names(rauldataset)
str(rauldataset)
tail(rauldataset)


all.equal(UMG, LDH3 + 24/LDH1, tolerance = 0.01)




##  - da detti valori di LDH1 e LDH3, calcolare l’indice di rischio UMG,
##  dove UMG = LDH3 + 24/LDH1;

LDH1inv = 24/LDH1
plot(LDH3, LDH1inv, col = "gray", pch =21)
points(LDH3[outcome == "malignant"], LDH1inv[outcome == "malignant"], col = "red")



library(ggplot2)
library(viridis)
library(dplyr)


# Partiamo dal dataset originale
df <- rauldataset %>%
  mutate(
    outcome = factor(outcome, levels = c("benign", "malignant")),
    # Calcolo del reciproco come da tuo esempio: LDH1inv = 24 / LDH1
    # Gestione sicura: se LDH1 <= 0 o NA, mettiamo NA in LDH1inv per evitare infiniti o errori.
    LDH1inv = ifelse(!is.na(LDH1) & LDH1 > 0, 24 / LDH1, NA_real_)
  )

# Opzionale: rimuoviamo righe con NA nelle variabili di interesse
df_plot <- df %>% filter(!is.na(LDH3), !is.na(LDH1inv), !is.na(outcome))

p <- ggplot(df_plot, aes(x = LDH3, y = LDH1inv,
                         color = outcome, shape = outcome)) +
  # punti come richiesto
  geom_point(alpha = 0.9, size = 2) +
  
  # palette viridis: benign = green bullet, malignant = purple triangle
  scale_color_viridis(
    discrete = TRUE,
    option = "D",
    begin  = 0.75,  # verde
    end    = 0.15   # viola
  ) +
  # benign bullet (16), malignant triangle (17)
  scale_shape_manual(values = c(16, 17)) +
  
  labs(
    x = "LDH3",
    y = expression("24 / LDH1"),
    color = "Outcome",
    shape = "Outcome",
    title = "(24 / LDH1) vs. LDH3",
    subtitle = " "
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p


## boxplot

library(ggplot2)
library(viridis)

df <- rauldataset
df$outcome <- factor(df$outcome, levels = c("benign", "malignant"))

p_box <- ggplot(df, aes(x = outcome, y = UMG, fill = outcome)) +
  
  # boxplot elegante
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.alpha = 0.6) +
  
  # palette coerente con gli altri plot: benign verde, malignant viola
  scale_fill_viridis(
    discrete = TRUE,
    option = "D",
    begin = 0.75,   # verde
    end   = 0.15    # viola
  ) +
  
  labs(
    x = "",
    y = "UMG index",
    title = "UMG index vs. outcome",
    subtitle = " ",
    fill = "Outcome"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",  # chiaro dalla x-label chi è chi
    axis.text.x = element_text(size = 13)
  )

p_box  ## fine boxplot



##### 1. total LDH
####

##  - considerare il valore di LDHTOT;
##  - definire un primo indicatore, dove detto primo indicatore mette in
##  correlazione detto UMG e detto LDHTOT, laddove detto campione
##  è classificato come segue sulla base di detto primo indicatore:
##  UMG > 40 – 0.05 * LDHTOT classe c)
##  UMG = 40 – 0.05 * LDHTOT classe b)
##  UMG < 40 – 0.05 * LDHTOT classe a)



library(ggplot2)
library(viridis)

df <- rauldataset
df$outcome <- factor(df$outcome, levels = c("benign", "malignant"))

p <- ggplot(df, aes(x = LDHtot, y = UMG,
                    color = outcome,
                    shape = outcome)) +
  
  geom_point(alpha = 0.8, size = 2) +
  
  #  logarithmic
  scale_x_log10() +
  
  # palette viridis: benign = green bullet, malignant = purple triangle
  scale_color_viridis(
    discrete = TRUE,
    option = "D",
    begin = 0.75,   # green
    end   = 0.15    # purple
  ) +
  
  # benign bullet, malignant triangle
  scale_shape_manual(values = c(16, 17)) +
  
  labs(
    x = "total LDH (log10)",
    y = "UMG index",
    color = "Outcome",
    shape = "Outcome",
    title = "UMG index vs. total LDH relation"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )


# p


# Aggiunta della linea di demarcazione (retta nei dati originali)
# Nota: la scala log10 viene applicata solo all'asse x; la funzione è valutata su x "reali".
p_line <- p +
  geom_function(
    aes(linetype = "first classifier"),
    fun = function(x) 40 - 0.05 * x,
    color = "black",
    size = 0.9
  ) +
  scale_linetype_manual(
    name = "",  # legenda pulita per la linea
    values = c("first classifier" = "solid")
  )

p_line



p_line + coord_cartesian(ylim = c(12, 35)) +
  annotate("text", x = 550, y = 17,
           label = "UMG = 40 - 0.05 × LDHtot",
           hjust = 0, size = 4.5)






table(outcome, class1)

# Install only if needed
# install.packages(c("pROC", "epiR"))

library(pROC)


# Assicurati che i livelli siano nell'ordine corretto
rauldataset$outcome <- factor(rauldataset$outcome, levels = c("benign", "malignant"))
rauldataset$class1  <- factor(rauldataset$class1,  levels = c("A", "C"))

# Confusion matrix (righe = outcome reale, colonne = test predetto)
cm <- table(rauldataset$outcome, rauldataset$class1)
cm
#            Pred
# Truth    A      C
# benign   TN     FP
# malignant FN     TP

TN <- cm["benign","A"]
FP <- cm["benign","C"]
FN <- cm["malignant","A"]
TP <- cm["malignant","C"]

N  <- sum(cm)

# Metriche puntuali
sensitivity <- TP / (TP + FN)                    # True Positive Rate
specificity <- TN / (TN + FP)                    # True Negative Rate
ppv         <- TP / (TP + FP)                    # Precision
npv         <- TN / (TN + FN)
accuracy    <- (TP + TN) / N
f1          <- if ((2*TP + FP + FN) > 0) 2*TP / (2*TP + FP + FN) else NA

# Intervalli di confidenza esatti (binomiali)
sens_ci <- binom.test(TP, TP + FN)$conf.int
spec_ci <- binom.test(TN, TN + FP)$conf.int
ppv_ci  <- binom.test(TP, TP + FP)$conf.int
npv_ci  <- binom.test(TN, TN + FN)$conf.int
acc_ci  <- binom.test(TP + TN, N)$conf.int

# Riepilogo ordinato
diag_summary <- data.frame(
  Metric = c("Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "F1 score"),
  Estimate = c(sensitivity, specificity, ppv, npv, accuracy, f1),
  CI_low = c(sens_ci[1], spec_ci[1], ppv_ci[1], npv_ci[1], acc_ci[1], NA),
  CI_high= c(sens_ci[2], spec_ci[2], ppv_ci[2], npv_ci[2], acc_ci[2], NA)
)
diag_summary







# Outcome binario come 0/1 per pROC: 1 = "malignant"
y <- ifelse(rauldataset$outcome == "malignant", 1, 0)

# Predittore continuo: più alto = maggiore rischio
score <- rauldataset$classifier1

# Costruzione ROC
roc_obj <- roc(response = y, predictor = score, direction = ">")  # ">" assume score più alto = più probabile positivo
roc_obj

# AUC e CI (DeLong)
auc_value <- auc(roc_obj)
auc_ci    <- ci.auc(roc_obj, method = "delong")

auc_value
auc_ci

# Best threshold secondo Youden's J
best_coords <- coords(roc_obj, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"))
best_coords

# Plot ROC con pROC
plot(roc_obj, col = "#3B82F6", lwd = 3, main = "ROC curve — classifier1 (UMG–LDHtot function)")
abline(a = 0, b = 1, lty = 2, col = "gray50")
legend("bottomright",
       legend = paste0("AUC = ", round(as.numeric(auc_value), 3),
                       " (95% CI ", round(as.numeric(auc_ci[1]), 3), "–",
                       round(as.numeric(auc_ci[3]), 3), ")"),
       bty = "n")




# Tabella per il test binario (class1)
binary_perf <- data.frame(
  Measure = c("Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "F1"),
  Estimate = c(sensitivity, specificity, ppv, npv, accuracy, f1),
  CI_95 = c(
    sprintf("[%.3f–%.3f]", sens_ci[1], sens_ci[2]),
    sprintf("[%.3f–%.3f]", spec_ci[1], spec_ci[2]),
    sprintf("[%.3f–%.3f]", ppv_ci[1], ppv_ci[2]),
    sprintf("[%.3f–%.3f]", npv_ci[1], npv_ci[2]),
    sprintf("[%.3f–%.3f]", acc_ci[1], acc_ci[2]),
    "" # F1 senza CI qui
  )
)

# Riga con AUC della ROC continua
auc_row <- data.frame(
  Measure = "AUC (continuous classifier1)",
  Estimate = as.numeric(auc_value),
  CI_95 = sprintf("[%.3f–%.3f]", as.numeric(auc_ci[1]), as.numeric(auc_ci[3]))
)

rbind(binary_perf, auc_row)





####  2



library(ggplot2)
library(viridis)

df <- rauldataset
df$outcome <- factor(df$outcome, levels = c("benign", "malignant"))

# dataset limits
xmin <- min(df$UMG, na.rm = TRUE)
xmax <- max(df$UMG, na.rm = TRUE)
ymin <- min(df$LDH5, na.rm = TRUE)
ymax <- max(df$LDH5, na.rm = TRUE)

# ---- GREEN REGION (UMG < 29) ----
rect_green <- data.frame(
  xmin = xmin,
  xmax = 29,
  ymin = ymin,
  ymax = ymax
)

# ---- GRAY REGION  ----
#  A = (29, 1.5)
#  B = (29, 13.7)
#  C = (32.5, 13.7)

A <- c(29, 1.5)
B <- c(29, 13.7)
C <- c(32.5, 13.7)

poly_gray <- data.frame(
  x = c(A[1], B[1], C[1]),
  y = c(A[2], B[2], C[2])
)


# ----   plot ----

p <- ggplot(df, aes(x = UMG, y = LDH5)) +
  
  # red region high risk
  annotate("rect",
           xmin = xmin, xmax = xmax,
           ymin = ymin, ymax = ymax,
           fill = "red", alpha = 0.08) +
  
  # green region low risk
  geom_rect(data = rect_green,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "green", alpha = 0.08) +
  
  # gray region higher risk
  geom_polygon(data = poly_gray,
               aes(x = x, y = y),
               inherit.aes = FALSE,
               fill = "grey60", alpha = 0.20) +
  
  
  # ---- CONTORNO DELLA REGIONE GRIGIA ----
geom_polygon(
  data = poly_gray,
  aes(x = x, y = y),
  inherit.aes = FALSE,
  fill = NA,
  color = "black",
  linewidth = 0.8
) +
  
  # ---- DOTS ----
geom_point(aes(color = outcome, shape = outcome),
           alpha = 0.9, size = 2) +
  
  # Palette viridis: benign=green, malignant=purple
  scale_color_viridis(
    discrete = TRUE,
    option = "D",
    begin  = 0.75,  # green
    end    = 0.15   # purple
  ) +
  
  scale_shape_manual(values = c(16, 17)) +
  
  labs(
    x = "UMG index",
    y = "LDH5",
    color = "Outcome",
    shape = "Outcome",
    title = "LDH5 vs. UMG: risk regions"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p





table(outcome, class2)

library(pROC)

# Assicuriamo fattori nell'ordine corretto
rauldataset$outcome <- factor(rauldataset$outcome, levels = c("benign", "malignant"))
rauldataset$class2   <- factor(rauldataset$class2,  levels = c("A", "B", "C"))

# -----------------------------
# 1) CREAZIONE TEST BINARIO
# -----------------------------
# A = negativo
# B,C = positivo
rauldataset$class2_bin <- ifelse(rauldataset$class2 == "A", "neg", "pos")
rauldataset$class2_bin <- factor(rauldataset$class2_bin, levels = c("neg", "pos"))

# -----------------------------
# 2) CONFUSION MATRIX
# -----------------------------
cm2 <- table(rauldataset$outcome, rauldataset$class2_bin)
cm2

# Estrazione TN, FP, FN, TP
TN <- cm2["benign","neg"]
FP <- cm2["benign","pos"]
FN <- cm2["malignant","neg"]
TP <- cm2["malignant","pos"]

N  <- sum(cm2)

# -----------------------------
# 3) METRICHE DIAGNOSTICHE
# -----------------------------
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
ppv         <- TP / (TP + FP)
npv         <- TN / (TN + FN)
accuracy    <- (TP + TN) / N
f1          <- 2 * TP / (2*TP + FP + FN)

# CI esatti (binomiali)
sens_ci <- binom.test(TP, TP + FN)$conf.int
spec_ci <- binom.test(TN, TN + FP)$conf.int
ppv_ci  <- binom.test(TP, TP + FP)$conf.int
npv_ci  <- binom.test(TN, TN + FN)$conf.int
acc_ci  <- binom.test(TP + TN, N)$conf.int

# Tabella di riepilogo
diag_summary2 <- data.frame(
  Metric = c("Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "F1 score"),
  Estimate = c(sensitivity, specificity, ppv, npv, accuracy, f1),
  CI_low = c(sens_ci[1], spec_ci[1], ppv_ci[1], npv_ci[1], acc_ci[1], NA),
  CI_high= c(sens_ci[2], spec_ci[2], ppv_ci[2], npv_ci[2], acc_ci[2], NA)
)

diag_summary2


     

# Outcome numerico per pROC: 1 = "malignant"
y <- ifelse(rauldataset$outcome == "malignant", 1, 0)

# Predittore continuo: p (probabilità di malignità di Classifier 3)
score <- rauldataset$p

# ROC
roc_p <- roc(response = y, predictor = score, direction = ">")
roc_p

# AUC + intervallo di confidenza (DeLong)
auc_p  <- auc(roc_p)
auc_ci <- ci.auc(roc_p, method = "delong")

auc_p
auc_ci

# Soglia di Youden (per completezza statistica)
best_p <- coords(roc_p, "best", best.method = "youden",
                 ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"))
best_p







######## 3














xbenign = UMG[outcome == "benign"]
ybenign = LDH5[outcome == "benign"]
length(xbenign)
length(ybenign)
mx = mean(xbenign)		# 23.07888
sx = sd(xbenign)			# 2.390284
my = mean(ybenign)		# 9.267797
sy = sd(ybenign)			# 1.906641
r = cor(xbenign, ybenign)	# -0.2069279
n=1000
# genero due N(0,1) indipendenti
n1=rnorm(n)
n2=rnorm(n)
# combino
n3=r*n1+sqrt(1-r^2)*n2
# allora
xx = mx+sx*n1
yy = my+sy*n3
# sono le simulazioni cercate 
par(mfrow=c(1,2))
plot(xbenign, ybenign, col = "chartreuse4")
plot(xx,yy, col = "chartreuse4")

summary(UMG)

xmalignant = UMG[outcome == "malignant"]
ymalignant = LDH5[outcome == "malignant"]



u1 = runif(n, 29, 40)
# u2 = runif(n)
u3 = runif(n, 1, 40)

par(mfrow=c(1,2))
plot(xmalignant, ymalignant, col = "coral3")
plot(u1,u3, col = "coral3")
par(mfrow=c(1,1))



xtest = 26
ytest = 12


x = c(xx, u1)
y = c(yy, u3)

length(x)

giudizio = factor(c(rep("verde", n), rep("rosso", n)))

# cbind(x,y,giudizio)
distanza = sqrt( (x-xtest)^2 + (y-ytest)^2 )
#boxplot(distanza)
raggio = quantile(distanza, 0.1)
vicini = which(distanza < raggio)
# giudizio[vicini]
table(giudizio[vicini])




xtest = 29
ytest = 5


x = c(xx, u1)
y = c(yy, u3)

length(x)

giudizio = factor(c(rep("verde", n), rep("rosso", n)))

# cbind(x,y,giudizio)
distanza = sqrt( (x-xtest)^2 + (y-ytest)^2 )
#boxplot(distanza)
raggio = quantile(distanza, 0.1)
vicini = which(distanza < raggio)
# giudizio[vicini]
table(giudizio[vicini])





####

xgren = xx[1:400]
ygren = yy[1:400]
xreed = u1[1:400]
yreed = u3[1:400]


theta = seq(0, 2*pi, .01)
xcirc = 26 + 4.78 * cos(theta)
ycirc = 13 + 4.78 * sin(theta)
plot(NA, xlim = c(20,32), ylim = c(5,20), xlab = "UMG", ylab = "LDH5")
lines(c(26, 26), c(0, 13), lty = 2, lwd = 3, col = "darkturquoise")
lines(c(0, 26), c(13, 13), lty = 2, lwd = 3, col = "darkturquoise")
text(xgren, ygren,  "low", col = "dodgerblue4")
text(xreed, yreed, "high", col = "deepskyblue4")
points(26, 13, lwd = 6, pch = 23, col = "darkcyan", bg = "chartreuse3")
lines(xcirc, ycirc, lty = 3, lwd = 3, col = "blue4")


###



table(outcome, class3)

library(pROC)

# -----------------------------
# 1) Fattori nell’ordine corretto
# -----------------------------
rauldataset$outcome <- factor(rauldataset$outcome, levels = c("benign", "malignant"))
rauldataset$class3   <- factor(rauldataset$class3,  levels = c("A", "B", "C"))

# -----------------------------
# 2) CREAZIONE DEL TEST BINARIO
# -----------------------------
# Classe A = negativo
# Classi B e C = positivo
rauldataset$class3_bin <- ifelse(rauldataset$class3 == "A", "neg", "pos")
rauldataset$class3_bin <- factor(rauldataset$class3_bin, levels = c("neg", "pos"))

# -----------------------------
# 3) CONFUSION MATRIX
# -----------------------------
cm3 <- table(rauldataset$outcome, rauldataset$class3_bin)
cm3

# TN, FP, FN, TP
TN <- cm3["benign","neg"]
FP <- cm3["benign","pos"]
FN <- cm3["malignant","neg"]
TP <- cm3["malignant","pos"]

N <- sum(cm3)

# -----------------------------
# 4) METRICHE DIAGNOSTICHE
# -----------------------------
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
ppv         <- TP / (TP + FP)
npv         <- TN / (TN + FN)
accuracy    <- (TP + TN) / N
f1          <- 2 * TP / (2*TP + FP + FN)

# -----------------------------
# 5) INTERVALLI DI CONFIDENZA
# -----------------------------
sens_ci <- binom.test(TP, TP + FN)$conf.int
spec_ci <- binom.test(TN, TN + FP)$conf.int
ppv_ci  <- binom.test(TP, TP + FP)$conf.int
npv_ci  <- binom.test(TN, TN + FN)$conf.int
acc_ci  <- binom.test(TP + TN, N)$conf.int

# -----------------------------
# 6) TABELLA DI RIEPILOGO
# -----------------------------
diag_summary3 <- data.frame(
  Metric = c("Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "F1 score"),
  Estimate = c(sensitivity, specificity, ppv, npv, accuracy, f1),
  CI_low = c(sens_ci[1], spec_ci[1], ppv_ci[1], npv_ci[1], acc_ci[1], NA),
  CI_high= c(sens_ci[2], spec_ci[2], ppv_ci[2], npv_ci[2], acc_ci[2], NA)
)

diag_summary3




# B. ROC and AUC using the probability p (from classifier 3)

library(pROC)
y <- ifelse(rauldataset$outcome == "malignant", 0, 1)
roc_p <- roc(y, rauldataset$p, direction = ">")

plot(roc_p, col = "darkred", lwd = 3, main = "ROC curve ")
auc(roc_p)
ci.auc(roc_p)


## sempre class 3

# Outcome numerico per pROC: 1 = "malignant"
y <- ifelse(rauldataset$outcome == "malignant", 0, 1)

# Predittore continuo: p (probabilità di malignità di Classifier 3)
score <- rauldataset$p

# ROC
roc_p <- roc(response = y, predictor = score, direction = ">")
roc_p

# AUC + intervallo di confidenza (DeLong)
auc_p  <- auc(roc_p)
auc_ci <- ci.auc(roc_p, method = "delong")

auc_p
auc_ci

# Soglia di Youden (per completezza statistica)
best_p <- coords(roc_p, "best", best.method = "youden",
                 ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"))
best_p


# Punto ROC del test binario corrente (A vs B+C)
fpr_bin <- 1 - specificity
tpr_bin <- sensitivity

# Plot ROC
plot(roc_p, col = "#7C3AED", lwd = 3,
     main = "ROC curve — Classifier 3 (p as continuous score)")
abline(a = 0, b = 1, lty = 2, col = "gray50")

legend("bottomright",
       legend = paste0("AUC = ", round(as.numeric(auc_p), 3),
                       " (95% CI ", round(as.numeric(auc_ci[1]), 3), "–",
                       round(as.numeric(auc_ci[3]), 3), ")"),
       bty = "n")

# Sovrapponi il punto del cut-off binario A vs (B+C)
points(fpr_bin, tpr_bin, pch = 19, col = "#DC2626", cex = 1.3)
text(fpr_bin, tpr_bin, labels = " binary cut (A vs B+C)", pos = 4, col = "#DC2626")


# Plot ROC
plot(roc_p, col = "darkred", lwd = 3,
     main = "ROC curve — Classifier 3 (p as continuous score)")
abline(a = 0, b = 1, lty = 2, col = "gray50")

legend("bottomright",
       legend = paste0("AUC = ", round(as.numeric(auc_p), 3),
                       " (95% CI ", round(as.numeric(auc_ci[1]), 3), "–",
                       round(as.numeric(auc_ci[3]), 3), ")"),
       bty = "n")



# Tabella binaria con CI formattati
binary_perf3 <- data.frame(
  Measure  = c("Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "F1"),
  Estimate = c(sensitivity, specificity, ppv, npv, accuracy, f1),
  CI_95    = c(
    sprintf("[%.3f–%.3f]", sens_ci[1], sens_ci[2]),
    sprintf("[%.3f–%.3f]", spec_ci[1], spec_ci[2]),
    sprintf("[%.3f–%.3f]", ppv_ci[1],  ppv_ci[2]),
    sprintf("[%.3f–%.3f]", npv_ci[1],  npv_ci[2]),
    sprintf("[%.3f–%.3f]", acc_ci[1],  acc_ci[2]),
    ""  # F1 senza CI
  )
)

# Riga AUC della ROC continua (p)
auc_row3 <- data.frame(
  Measure  = "AUC (continuous p)",
  Estimate = as.numeric(auc_p),
  CI_95    = sprintf("[%.3f–%.3f]", as.numeric(auc_ci[1]), as.numeric(auc_ci[3]))
)

perf_table3 <- rbind(binary_perf3, auc_row3)
perf_table3








#### fine 
# "raul"    

table(outcome, raul)


rauldataset$raul_bin <- ifelse(rauldataset$raul == "A", "neg", "pos")
rauldataset$raul_bin <- factor(rauldataset$raul_bin, levels = c("neg", "pos"))

cm <- table(rauldataset$outcome, rauldataset$raul_bin)
cm

TN <- cm["benign","neg"]
FP <- cm["benign","pos"]
FN <- cm["malignant","neg"]
TP <- cm["malignant","pos"]

sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
ppv <- TP / (TP + FP)
npv <- TN / (TN + FN)
accuracy <- (TP + TN) / sum(cm)

cbind(sensitivity, specificity, ppv, npv, accuracy)





# C. Calibration analysis (optional but very impressive)

library(dplyr)

cal <- rauldataset %>%
  mutate(bin = cut(p, breaks = seq(0,1,length=11))) %>%
  group_by(bin) %>%
  summarise(mean_p = mean(p),
            obs_rate = mean(outcome == "malignant"),
            n = n())
cal





#  D. Concordance between classifiers: agreement analysis
# You can show that:
  
#   classifier 1 and classifier 2 agree strongly
# classifier 3 reinforces high‑risk areas
# RAUL consolidates their strengths
#  R code (Cohen’s kappa)



library(psych)
cohen.kappa(cbind(as.numeric(rauldataset$class1),
                  as.numeric(rauldataset$class2)))

cohen.kappa(cbind(as.numeric(rauldataset$class1),
                  as.numeric(rauldataset$class3)))


LR_pos <- sensitivity / (1 - specificity)
LR_neg <- (1 - sensitivity) / specificity

cbind(LR_pos, LR_neg)



library(ggplot2)

ggplot(rauldataset, aes(x = UMG, fill = raul)) +
  geom_density(alpha = 0.6)




##%%##%%##%%
n=1000




xbenign = UMG[outcome == "benign"]
ybenign = LDH5[outcome == "benign"]
mx = mean(xbenign); sx = sd(xbenign); 
my = mean(ybenign); sy = sd(ybenign);
r = cor(xbenign, ybenign)


n1=rnorm(n); n2=rnorm(n)
n3=r*n1+sqrt(1-r^2)*n2


xB = mx+sx*n1
yB = my+sy*n3

xM = runif(n, 29, 40)
yM = runif(n, 1, 40)

x = c(xB, xM)
y = c(yB, yM)
colour = factor(c(rep("green", n), rep("purple", n)))
plot(x,y, col = colour)




xtest = 28.5
ytest = 13.5
colour = factor(c(rep("green", n), rep("purple", n)))
distance = sqrt( (x-xtest)^2 + (y-ytest)^2 )
radius = quantile(distance, 0.1)
nearest = which(distance < radius)
table(colour[nearest])
(p = table(colour[nearest])[[2]] / sum(table(colour[nearest])))



library(ggplot2)
library(viridis)

set.seed(123)  # per riproducibilità (opzionale)
n <- 1000

# --- Dati di partenza: estraiamo le statistiche sui benign come nel tuo script ---
xbenign <- rauldataset$UMG[rauldataset$outcome == "benign"]
ybenign <- rauldataset$LDH5[rauldataset$outcome == "benign"]

mx <- mean(xbenign, na.rm = TRUE); sx <- sd(xbenign, na.rm = TRUE)
my <- mean(ybenign, na.rm = TRUE); sy <- sd(ybenign, na.rm = TRUE)
r  <- cor(xbenign, ybenign, use = "complete.obs")

# --- Simulazione bivariata con correlazione r per il cluster "green" ---
n1 <- rnorm(n); n2 <- rnorm(n)
n3 <- r * n1 + sqrt(1 - r^2) * n2

xB <- mx + sx * n1
yB <- my + sy * n3

# --- Cluster "purple" uniforme nel rettangolo specificato ---
xM <- runif(n, 29, 40)
yM <- runif(n, 1, 40)

# --- Costruzione del data.frame per ggplot ---
x <- c(xB, xM)
y <- c(yB, yM)
colour <- factor(c(rep("green", n), rep("purple", n)), levels = c("green", "purple"))

df_plot <- data.frame(
  x = x,
  y = y,
  group = colour
)

# --- Plot ggplot2 (coerente con lo stile precedente) ---
p <- ggplot(df_plot, aes(x = x, y = y, color = group, shape = group)) +
  geom_point(alpha = 0.85, size = 1.8) +
  
  # Palette viridis: primo livello (green) -> verde, secondo (purple) -> viola
  scale_color_viridis(
    discrete = TRUE,
    option = "D",
    begin  = 0.75,  # verde
    end    = 0.15   # viola
  ) +
  # forme: green = cerchio, purple = triangolo
  scale_shape_manual(values = c(16, 17)) +
  
  labs(
    x = "x (simulated UMG)",
    y = "y (simulated LDH5)",
    color = "Group",
    shape = "Group",
    title = "Simulated clusters (green vs purple)",
    subtitle = "benign-like cluster (gaussian); malignant-like cluster (uniform)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

p

p +
  stat_density_2d(data = subset(df_plot, group == "green"),
                  aes(fill = after_stat(level)),
                  geom = "polygon", bins = 5,
                  alpha = 0.12, color = NA, show.legend = FALSE)


##$$$%%%$$$%$%$%$


library(ggplot2)
library(viridis)

# df_plot già definito (x, y, group)
df_plot$group <- factor(df_plot$group, levels = c("green", "purple"))

# Punto di test e distanze
xtest <- 28.5
ytest <- 13.5
df_plot$distance <- sqrt((df_plot$x - xtest)^2 + (df_plot$y - ytest)^2)

# Raggio = 10° percentile (≈200 su 2000 punti totali)
k <- ceiling(0.10 * nrow(df_plot))
radius <- sort(df_plot$distance)[k]

# Selezione dei "nearest"
df_plot$nearest <- df_plot$distance < radius  # usa <= se vuoi includere i tie
df_near <- subset(df_plot, nearest)

# Circonferenza
theta <- seq(0, 2*pi, length.out = 361)
circle_df <- data.frame(
  x = xtest + radius * cos(theta),
  y = ytest + radius * sin(theta)
)

# Plot zoom centrato su (xtest, ytest)
p_zoom <- ggplot() +
  # (facoltativo) contesto dei punti totali in grigio chiaro:
   geom_point(data = df_plot, aes(x = x, y = y), color = "grey85", size = 1.2, alpha = 0.35) +
  
  geom_point(data = df_near,
             aes(x = x, y = y, color = group, shape = group),
             alpha = 0.9, size = 2) +
  
  geom_path(data = circle_df, aes(x = x, y = y),
            color = "black", linewidth = 0.7) +
  
  annotate("point", x = xtest, y = ytest, shape = 4, size = 3.5,
           stroke = 1.1, colour = "black") +
  
  coord_cartesian(xlim = c(xtest - radius, xtest + radius),
                  ylim = c(ytest - radius, ytest + radius),
                  expand = FALSE) +
  
  scale_color_viridis(discrete = TRUE, option = "D", begin = 0.75, end = 0.15) +
  scale_shape_manual(values = c(16, 17)) +
  
  labs(
    x = "x (simulated UMG)",
    y = "y (simulated LDH5)",
    color = "Group",
    shape = "Group",
    title = "Zoom around test point"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

p_zoom

p_zoom <- p_zoom + theme(legend.position = "none")
