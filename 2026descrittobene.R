

set.seed(123)
www = "https://raw.githubusercontent.com/umgstat/raul/refs/heads/main/rauldataset.csv"
rauldataset = read.csv(www, header = TRUE, stringsAsFactors = TRUE)
attach(rauldataset)
head(rauldataset)

names(rauldataset)
str(rauldataset)
tail(rauldataset)






##### 1. total LDH

##  - da detti valori di LDH1 e LDH3, calcolare l’indice di rischio UMG,
##  dove UMG = LDH3 + 24/LDH1;
##  - considerare il valore di LDHTOT;
##  - definire un primo indicatore, dove detto primo indicatore mette in
##  correlazione detto UMG e detto LDHTOT, laddove detto campione
##  è classificato come segue sulla base di detto primo indicatore:
##  UMG > 40 – 0.05 * LDHTOT classe c)
##  UMG = 40 – 0.05 * LDHTOT classe b)
##  UMG < 40 – 0.05 * LDHTOT classe a)


min(UMG - (LDH3 + 24/LDH1))
max(UMG - (LDH3 + 24/LDH1))



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
  annotate("text", x = 700, y = 17,
           label = "UMG = 40 - 0.05 × LDHtot",
           hjust = 0, size = 4.5)






####



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







########














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





xsane = UMG[OUTCOME == "benignant"]
ysane= LDH5[OUTCOME == "benignant"]
length(xsane)
length(ysane)
mx = mean(xsane)		# 23.07888
sx = sd(xsane)			# 2.390284
my = mean(ysane)		# 9.267797
sy = sd(ysane)			# 1.906641
r = cor(xsane, ysane)	# -0.2069279
n=50000
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
plot(xsane, ysane, col = "chartreuse4")
plot(xx,yy, col = "chartreuse4")

summary(UMG)

xmalate = UMG[OUTCOME == "malignant"]
ymalate = LDH5[OUTCOME == "malignant"]

summary(xmalate)
summary(ymalate)


u1 = runif(n, 29, 40)
# u2 = runif(n)
u3 = runif(n, 1, 40)

par(mfrow=c(1,2))
plot(xmalate, ymalate, col = "coral3")
plot(u1,u3, col = "coral3")




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


