
# Fit log acute costs to log healthcare expenditure per capita

library("dplyr")
library("data.table")
library("ggplot2")

# read data
data.costs <- read.csv("./data/data.acute_costs_2020USD.csv")
data.countries <- read.csv("./data/data.countries.csv")
dataset <- as.data.table(inner_join(data.costs,data.countries,by="iso3"))

dataset <- dataset[study_year>2004] # remove older studies from USA
dataset[, log_wb_health_expenditure_per_capita_usd_2019 := log(wb_health_expenditure_per_capita_usd_2019)]

# linear fit
lm.3 <-  lm(log(average_cost_2020_USD) ~ log_wb_health_expenditure_per_capita_usd_2019, data=dataset)

#predict per country
results <- data.table(
  country = data.countries$iso3,
  wb_gdp_per_capita_usd_2020 = data.countries$wb_gdp_per_capita_usd_2020,
  wb_health_expenditure_per_capita_usd_2019 = data.countries$wb_health_expenditure_per_capita_usd_2019
)[order(wb_health_expenditure_per_capita_usd_2019)]

results[, log_wb_health_expenditure_per_capita_usd_2019 := log(wb_health_expenditure_per_capita_usd_2019)]

results <- cbind(
  results[,.(
    country,
    wb_health_expenditure_per_capita_usd_2019,
    log_wb_health_expenditure_per_capita_usd_2019,
    wb_gdp_per_capita_usd_2020
  )],
  exp(predict(lm.3,results,interval="confidence"))
)

# plot <- ggplot(data=results, aes(x = wb_health_expenditure_per_capita_usd_2019)) +
#   geom_line(aes(y=fit), color="blue") +
#   geom_ribbon(aes(ymin=lwr, ymax=upr),alpha=0.4) +
#   geom_point(data = dataset[,country:=iso3],
#              aes(x=wb_health_expenditure_per_capita_usd_2019, y=average_cost_2020_USD),
#              color="red", 
#              shape="square"
#   ) +
#   ylab("Cost of acute care (2020 USD)") +
#   xlab("Per capita healthcare expenditure (2019 USD)") +
#   ggtitle("Regression of log(costs) against log(per capita healthcare exenditure)") +
#   NULL
# 
# plot.log <- plot + scale_y_log10()
# plot.log_log <- plot + scale_y_log10() + scale_x_log10()

# export result
pred <- predict(lm.3,results,se.fit = TRUE)
results.out <- cbind(country = results$country, fit = pred$fit, se = pred$se.fit)
# fwrite(results.out,"acute_cost_fit.csv")



 

