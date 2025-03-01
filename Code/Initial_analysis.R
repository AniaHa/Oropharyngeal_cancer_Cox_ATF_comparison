# paczki
library(survminer)
library(foreign)
library(survival)
library(KMsurv)
library(survMisc) 

# wczytanie danych 

dane <- read.csv("E:/Studia/PW/II SEMESTR/Biostatystyka/projekty/Proj 5/pharynx.csv")
head(dane)

dane$INST <- as.factor(dane$INST)
dane$SEX <- as.factor(dane$SEX)
dane$TX <- as.factor(dane$TX)
dane$GRADE <- as.factor(dane$GRADE)
dane$COND <- as.factor(dane$COND)
dane$SITE <- as.factor(dane$SITE)
dane$T_STAGE <- as.factor(dane$T_STAGE)
dane$N_STAGE <- as.factor(dane$N_STAGE)
dane$STATUS <- as.factor(dane$STATUS)
# age - ciagla
summary(dane)


# Podpunkt 1
# sprawdzamy wplyw zmiennych: SEX, T_STAGE, N_STAGE

# Funkcja przezycia pacjentow, bez wplywu zmiennych

result <- survfit(Surv(TIME, STATUS)~1, data = dane)
surv_median(result)
# fit1.png
ggsurvplot(
  result,
  data = dane,
  xlab = "Days",
  ylab = "Overall survival probability")

Estimated_prob<-cbind(result[["time"]],result[["surv"]])

###### SEX ##### 

# wp³yw sex na czas prze¿ycia
result.s <- survfit(Surv(TIME, STATUS)~SEX, data = dane)
# fit1.png
ggsurvplot(
  result.s,
  data = dane,
  xlab = "Days",
  ylab = "Overall survival probability",
  legend.title = "Sex",
  conf.int = TRUE,
  legend.labs = c("Man", "Woman"))

# krzywe prze¿ycia s¹ po³o¿one bardzo blisko siebie, 
# a przedzia³y ufnoœci dla poszczególnych warstw pokrywaj¹ siê
# w du¿ym stopniu. Mo¿emy zatem podejrzewaæ, 
# ¿e plec nie ma istotnego wp³ywu na szansê prze¿ycia

survdiff(Surv(TIME, STATUS)~SEX, data = dane)

# du¿e p-value - nie ma istotnych ró¿nic miêdzy p³ciami w przezyciu
survdiff(Surv(TIME, STATUS)~T_STAGE, data = dane)
t1 <- ten(Surv(TIME, STATUS)~SEX, data = dane)
comp(t1)
attributes(t1)$sup[1,]

###### STAGE #####   
 
# wplyw T_STAGE
dane3 <- dane2
dane3$T_STAGE[which(dane2$T_STAGE == 1)] <- 2

result.t <- survfit(Surv(TIME, STATUS)~T_STAGE, data = dane)
# fit1.png
ggsurvplot(
  result.t,
  data = dane,
  xlab = "Days",
  ylab = "Overall survival probability",
  legend.title = "T_STAGE",
  conf.int = TRUE,
  legend.labs = c("1", "2", "3", "4"))

survdiff(Surv(TIME, STATUS)~T_STAGE, data = dane)
t1 <- ten(Surv(TIME, STATUS)~T_STAGE, data = dane)
comp(t1)
attributes(t1)$tft$tft[1,]$Z

t_stage.test <- survdiff(Surv(TIME, STATUS)~T_STAGE, data = dane3)
l <- sum((t_stage.test$obs - t_stage.test$exp)*c(3,2,1))
m <- sqrt(sum(outer(1:3, 1:3, "*")*t_stage.test$var))
Z <- l/m # test statistic for testing trend, based on log-rank 
Z

##### N_STAGE #####
result.n <- survfit(Surv(TIME, STATUS)~N_STAGE, data = dane)
# fit1.png
ggsurvplot(
  result.n,
  data = dane,
  xlab = "Days",
  ylab = "Overall survival probability",
  legend.title = "N_STAGE",
  conf.int = TRUE,
  legend.labs = c("1", "2", "3", "4"))
# grupy 1 i 4 oraz 2 i 3 siê zlewaj¹

survdiff(Surv(TIME, STATUS)~N_STAGE, data = dane)


ggcoxdiagnostics(result.step, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())

cox.fit <- coxph(Surv(TIME, STATUS)~ COND, data = dane5)
plot(survfit(cox.fit,stype=1),mark="+")
plot(survfit(cox.fit, newdata = data.frame(COND=as.factor(c(1,2, 3)))), lty=c(1,2) , fun="cloglog", main=expression(log(-log(S(t)))))
legend("topright", legend=c("age = 50", "age = 60", "COND = 3"), lty=c(1,2, 3))

cox.fit <- coxph(Surv(TIME, STATUS)~ SEX, data = dane5)
plot(survfit(cox.fit, newdata = data.frame(SEX=as.factor(c(1,2)))), lty=c(1,2) , fun="cloglog", main=expression(log(-log(S(t)))))
legend("topright", legend=c("age = 50", "age = 60", "COND = 3"), lty=c(1,2, 3))

cox.fit <- coxph(Surv(TIME, STATUS)~ T_STAGE, data = dane5)
plot(survfit(cox.fit, newdata = data.frame(T_STAGE=as.factor(c(2, 3, 4)))), lty=c(1,2) , fun="cloglog", main=expression(log(-log(S(t)))))
legend("topright", legend=c("age = 50", "age = 60", "COND = 3"), lty=c(1,2, 3))

cox.fit <- coxph(Surv(TIME, STATUS)~ T_STAGE + SEX + COND, data = dane5)
plot(survfit(cox.fit, newdata = data.frame(T_STAGE=as.factor(c(2, 2, 2)), SEX = as.factor(c(1, 1, 1)),
             COND = as.factor(c(1, 2, 3)))), lty=c(1,2, 3) , fun="cloglog", main=expression(log(-log(S(t)))))
legend("topright", legend=c("age = 50", "age = 60", "COND = 3"), lty=c(1,2, 3))
plot(cox.fit)

#####################################
# model Coxa dla danych binarnych 

dane0 <- dane

#===== SITE =====
dane0$site_1 <- as.numeric(dane0$SITE == 1)
dane0$site_2 <- as.numeric(dane0$SITE == 2)

dane0 <- subset(dane0, select = -c(SITE))

#===== GRADE =====
dane0$grade_1 <- as.numeric(dane$grade == 1)
dane0$grade_2 <- as.numeric(dane$grade == 2)

dane0 <- subset(dane0, select = -c(grade))

#===== COND =====
dane0$cond_1 <- as.numeric(dane$cond == 1)
dane0$cond_2 <- as.numeric(dane$cond == 2)

dane0 <- subset(dane0, select = -c(cond))

#===== SEX =====
dane0$sex <- dane0$sex - 1

#===== TX =====
#0 - radioterapia
#1 - radio i chemo

dane0$tx <- dane0$tx - 1

