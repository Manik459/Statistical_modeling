
library(BayesLCA)
library(jsonlite)
library(mclust)
library(zoo)
library(MCMCpack)
library(readr)
library("ggplot2")
library("readr")
library(reshape)
library(plotly)

df_bus_TO <- stream_in(url("https://www.scss.tcd.ie/~arwhite/Teaching/CS7DS3/business_open_Toronto.json"))

df_bus_TO <- flatten(df_bus_TO)
head(df_bus_TO)
df_bus_TO$neighborhood <- factor(df_bus_TO$neighborhood)
nlevels(df_bus_TO$neighborhood)

df_bus_TO$neighborhood <- as.numeric(df_bus_TO$neighborhood)
df_bus_TO = df_bus_TO[!(df_bus_TO$neighborhood ==38),]
df_bus_TO = df_bus_TO[!(df_bus_TO$neighborhood ==14),]
agg <- aggregate(df_bus_TO[,10],list(df_bus_TO$neighborhood),mean)

ggplot(df_bus_TO) + geom_boxplot(aes(x = reorder(neighborhood, stars, median), stars, fill = reorder(neighborhood, stars, median)), show.legend=FALSE)
#stars Vs count 
ggplot(df_bus_TO, aes(stars)) + stat_bin()
# neighbourhood vs count
ggplot(df_bus_TO, aes(x = reorder(neighborhood, neighborhood, length))) + stat_count()
compare_m_gibbs <- function(y, ind, maxiter = 5000)
{
  
  ### weakly informative priors
  a0 <- 1.9 ; b0 <- 1 ## tau_w hyperparameters
  eta0 <-1/2 ; t0 <- 5 ## tau_b hyperparameters
  mu0<-3.5 ; gamma0 <- 1.25
  ###
  
  ### starting values
  m <- nlevels(ind)
  ybar <- theta <- tapply(y, ind, mean)
  tau_w <- mean(1 / tapply(y, ind, var),na.rm = TRUE) ##within group precision
  mu <- mean(theta)
  tau_b <-var(theta) ##between group precision
  n_m <- tapply(y, ind, length)
  an <- a0 + sum(n_m)/2
  ###
  
  ### setup MCMC
  theta_mat <- matrix(0, nrow=maxiter, ncol=m)
  mat_store <- matrix(0, nrow=maxiter, ncol=3)
  ###
  
  ### MCMC algorithm
  for(s in 1:maxiter) 
  {
    
    # sample new values of the thetas
    for(j in 1:m) 
    {
      taun <- n_m[j] * tau_w + tau_b
      thetan <- (ybar[j] * n_m[j] * tau_w + mu * tau_b) / taun
      theta[j]<-rnorm(1, thetan, 1/sqrt(taun))
    }
    
    #sample new value of tau_w
    ss <- 0
    for(j in 1:m){
      ss <- ss + sum((y[ind == j] - theta[j])^2)
    }
    bn <- b0 + ss/2
    tau_w <- rgamma(1, an, bn)
    
    #sample a new value of mu
    gammam <- m * tau_b + gamma0
    mum <- (mean(theta) * m * tau_b + mu0 * gamma0) / gammam
    mu <- rnorm(1, mum, 1/ sqrt(gammam)) 
    
    # sample a new value of tau_b
    etam <- eta0 + m/2
    tm <- t0 + sum((theta-mu)^2)/2
    tau_b <- rgamma(1, etam, tm)
    
    #store results
    theta_mat[s,] <- theta
    mat_store[s, ] <- c(mu, tau_w, tau_b)
  }
  colnames(mat_store) <- c("mu", "tau_w", "tau_b")
  return(list(params = mat_store, theta = theta_mat))
}
fit2 <- compare_m_gibbs(df_bus_TO$stars, factor(df_bus_TO$neighborhood))
apply(fit2$params, 2, mean)
apply(fit2$params, 2, sd)
mean(1/sqrt(fit2$params[, 3]))
sd(1/sqrt(fit2$params[, 3]))
#plot of the the clusters
datasim=as.data.frame(fit2$theta)
datasim

datasim$ID <- seq.int(nrow(datasim))


md <- melt(datasim,id="ID")

ggplot(md) + 
  geom_boxplot(aes(x = reorder(variable, value, mean), value, fill = reorder(variable, value, mean)), show.legend=TRUE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# comparing two neighbourhoods
ggplot(df_bus_TO) + geom_boxplot(aes(neighborhood, stars, fill = neighborhood)) + geom_jitter(aes(neighborhood, stars, shape = df_bus_TO$neighborhood))



sel_nbhd <- df_bus_TO$neighborhood
df_cht_kt <- df_bus_TO[sel_nbhd, ]


class(df_cht_kt$categories)

#  try to make tidy -inelegant solution here
bus_cat_tidy <- cbind(df_cht_kt[1, ]$business_id, df_cht_kt[1, ]$categories[[1]])
for(i in 2:nrow(df_cht_kt)) bus_cat_tidy <- rbind(bus_cat_tidy, cbind(df_cht_kt[i, ]$business_id, df_cht_kt[i, ]$categories[[1]]))
cat_names <- names(sort(table(bus_cat_tidy[, 2]), decreasing = TRUE))[2:10]

cat_bus_ind_mat <- t(sapply(tapply(bus_cat_tidy[, 2], bus_cat_tidy[, 1], function(y) cat_names %in% y), function(x) x*1))

colnames(cat_bus_ind_mat) <- cat_names
df_cat <- data.frame(ind = rownames(cat_bus_ind_mat), cat_bus_ind_mat)
business_merge <- merge(df_cht_kt, df_cat, by.x = "business_id", by.y = "ind")


##question 2

data1 <- read_rds("Final.rds")
fit_mcmc <- MCMCregress(data1$stars.y ~data1$useful+data1$funny+data1$cool+data1$review_count+data1$attributes.RestaurantsPriceRange2+data1$Food+data1$Nightlife+data1$Bars+data1$Sandwiches+data1$Breakfast...Brunch+data1$Chinese+data1$Canadian..New.+data1$Cafes+data1$Coffee...Tea)
dummy <- data1
dummy <- dummy[ -c(1,2:5) ]
dummy
drops=c("review_count" ,
        "useful" ,
        "cool" ,
        "PriceRange" ,
        "Canadian..New." ,
        "Chinese" ,
        "Coffee...Tea" ,
        "Sandwiches" ,
        "Nightlife" ,
        "Bars" ,
        "Cafes" ,
        "Breakfast...Brunch" ,
        "Food" ,
        "funny",
        "stars.y" )
dummy=dummy[,(names(dummy) %in% drops)]

fit_mcmc
par(mar=c(2,2,2,2))
summary(fit_mcmc)
plot(fit_mcmc)
data1$attributes.RestaurantsPriceRange2 = round(data1$attributes.RestaurantsPriceRange2)
beta_mean <- apply(fit_mcmc, 2, mean)
coeff <-as.data.frame(beta_mean)
coeff[-1] <- coeff[order(coeff$beta_mean),-1]
coeff
data1$attributes.RestaurantsPriceRange2 <- as.numeric(data1$attributes.RestaurantsPriceRange2)

df_dummy <- dummy[ ,-4]  # remove response variable
#df_dummy$sex <- as.numeric(df_dummy$sex) # for convenience, change class of sex variable - trickier if variable has several levels
strdf_dummy <- cbind(1, as.matrix(df_dummy)) # add dummy variable for intercept, convert to matrix class
foo <-beta_mean[-16]
pred_fit <- strdf_dummy %*% as.matrix(foo) # make prediction, ignore variance parameter
p <-plot(data1$stars.y,pred_fit ) 
data1$predicted <-pred_fit
p + geom_jitter()
plot_external <- ggplot(data1)+
  aes(x = data1$stars.y)+
  # fitted values
  geom_jitter(aes(y=predicted), color = "purple")+
  geom_point(aes(y=stars.y), color = "red")
plot_external
 
#aic bic check
lm1 <- lm(data1$stars.y ~(data1$useful+data1$funny+data1$cool+data1$review_count+data1$Food+data1$Nightlife+data1$Bars+data1$Sandwiches+data1$Breakfast...Brunch+data1$Chinese+data1$Canadian..New.+data1$Cafes+data1$Coffee...Tea)^2)
summary(lm1)
step_AIC2 <- step(lm(y ~ (sex + bmi + map + tc + ldl + ltg)^2, data = diabetes_df))
aic <- step(lm1)

step_BIC <- step(lm1, k=log(nrow(data1)))
dv <- dv[ -c(5,5:12) ]
dv$stars.x<- NULL
lm2 <- lm(dv$stars.y~.,dv)
aic <- step(lm2)
step_BIC <- step(lm2, k=log(nrow(dv)))

# fitting the McMc for the best variables
fit2 <- MCMCregress(data1$stars.y ~data1$useful+data1$funny+data1$cool+data1$review_count+data1$Food+data1$Nightlife+data1$Bars+data1$Sandwiches+data1$Breakfast...Brunch+data1$Chinese+data1$Canadian..New.+data1$Cafes+data1$Coffee...Tea, data = data1, B0 = 1, marginal.likelihood = "Chib95")
attr(fit2, "logmarglike")

fit3 <- MCMCregress(data1$stars.y ~data1$useful+data1$funny+data1$Food+data1$Nightlife+data1$Bars+data1$Sandwiches+data1$Breakfast...Brunch+data1$Chinese+data1$Canadian..New.+data1$Cafes+data1$Coffee...Tea, data = data1, B0 = 1, marginal.likelihood = "Chib95")
attr(fit3,"logmarglike")

##question 3
df_cht_kt=fromJSON("https://www.scss.tcd.ie/~arwhite/Teaching/CS7DS3/business_open_Toronto.json")


class(df_cht_kt$categories)
bus_cat_tidy <- cbind(df_cht_kt[1, ]$business_id, df_cht_kt[1, ]$categories[[1]])
for(i in 2:nrow(df_cht_kt)) bus_cat_tidy <- rbind(bus_cat_tidy, cbind(df_cht_kt[i, ]$business_id, df_cht_kt[i, ]$categories[[1]]))
cat_names <- names(sort(table(bus_cat_tidy[, 2]), decreasing = TRUE))[2:74]
cat_bus_ind_mat <- t(sapply(tapply(bus_cat_tidy[, 2], bus_cat_tidy[, 1], function(y) cat_names %in% y), function(x) x*1))


colnames(cat_bus_ind_mat) <- cat_names
df_cat <- data.frame(ind = rownames(cat_bus_ind_mat), cat_bus_ind_mat)
business_merge <- merge(df_cht_kt, df_cat, by.x = "business_id", by.y = "ind")


business_merge$attributes.RestaurantsPriceRange2 = na.aggregate(business_merge$attributes.RestaurantsPriceRange2 )
sum(is.na(business_merge))

drops=c("categories")
df_cat=business_merge[,!(names(business_merge) %in% drops)]


colnames(df_cat)[colnames(df_cat)=="attributes.RestaurantsPriceRange2"] <- "PriceRange"


drops=c("business_id",
        "review_id",
        "user_id",
        "stars.x",
        "date",
        "text",
        "useful",
        "funny",
        "cool",
        "address",
        "city",
        "name",
        "state",
        "postal_code",
        "latitude",
        "longitude",
        "stars",
        "review_count",
        "is_open",
        "PriceRange"
)
df2=df_cat[,!(names(df_cat) %in% drops)]

apply(df2[,-1], 2, mean)

length(df2)



df3=df2[1:(length(df2)-61)]


keep=c("Downtown Core",
       "Scarborough",
       "Etobicoke",
       "Entertainment District",
       "Willowdale",
       "The Danforth",
       "Leslieville",
       "Financial District",
       "Milliken",
       "Yorkville",
       "Mount Pleasant and Davisville",
       "St. Lawrence",
       "Chinatown",
       "Church-Wellesley Village",
       "Kensington Market",
       "The Annex",
       "The Junction",
       "Yonge and Eglinton",
       "Corktown",
       "Riverdale"
)

df4=df3[df3$neighborhood %in% keep ,]


fit_4 <- blca.em(df4[-1],4)
fit_5 <- blca.em(df4[-1],5)
fit_8 <- blca.em(df4[-1],8)
fit_6 <- blca.em(df4[-1],6)

fit_8$BIC
fit_5$BIC
fit_4$BIC
fit_6$BIC
plot(fit_5, which = 1)
plot(fit_5, which = 2)

plot(fit_5, which = 5)


table(MAP(Zscore(df4[-1], fit_5)), df4[,1])
table(MAP(Zscore(df4[-1], fit_5)))

plot(table(MAP(Zscore(df4[-1], fit_5)), df4[,1]), las=1)
plot(table(MAP(Zscore(df4[-1], fit_6)), df4[,1]), las=1)
