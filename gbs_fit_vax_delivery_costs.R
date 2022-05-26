
# Fit log delivery costs to log gdp per capita

library("dplyr")
library("data.table")
library("ggplot2")

# read data
data.costs <- read.csv("./data/data.vax_delivery_costs_2020USD.csv")
data.countries <- read.csv("./data/data.countries.csv")
dataset <- as.data.table(inner_join(data.costs,data.countries,by="iso3"))

dataset[, log_wb_gdp_per_capita_usd_2020 := log(wb_gdp_per_capita_usd_2020)]

# linear fit
lm.3 <-  lm(log(cost_per_dose_usd_2020) ~ log_wb_gdp_per_capita_usd_2020, data=dataset)

#predict per country
results <- data.table(
  country = data.countries$iso3,
  wb_gdp_per_capita_usd_2020 = data.countries$wb_gdp_per_capita_usd_2020,
  wb_health_expenditure_per_capita_usd_2019 = data.countries$wb_health_expenditure_per_capita_usd_2019
)[order(wb_gdp_per_capita_usd_2020)]

results[, log_wb_gdp_per_capita_usd_2020 := log(wb_gdp_per_capita_usd_2020)]

results <- cbind(
  results[,.(
    country,
    wb_health_expenditure_per_capita_usd_2019,
    log_wb_gdp_per_capita_usd_2020,
    wb_gdp_per_capita_usd_2020
  )],
  exp(predict(lm.3,results,interval="confidence"))
)

# plot <- ggplot(data=results, aes(x = wb_gdp_per_capita_usd_2020)) +
#   geom_line(aes(y=fit), color="blue") +
#   geom_ribbon(aes(ymin=lwr, ymax=upr),alpha=0.4) +
#   geom_point(data = dataset[,country:=iso3],
#              aes(x=wb_gdp_per_capita_usd_2020, y=cost_per_dose_usd_2020),
#              color="red", 
#              shape="square"
#   ) +
#   
#   ylab("Cost of vaccine delivery (2020 USD)") +
#   xlab("GDP per capita (2020 USD)") +
#   ggtitle("Regression of log(vax delivery costs) against log(GDP per capita)") +
#   NULL
# 
# plot.log <- plot + scale_y_log10()
# plot.log_log <- plot + scale_y_log10() + scale_x_log10()

# export result
pred <- predict(lm.3,results,se.fit = TRUE)
results.out <- cbind(country = results$country, fit = pred$fit, se = pred$se.fit)
# fwrite(results.out,"vax_delivery_cost_fit.csv")



 

