## Custom functions
fct_case_when <- function(...) {
  # uses dplyr::case_when but converts the output to a factor,
  # with factors ordered as they appear in the case_when's  ... argument
  args <- as.list(match.call())
  levels <- sapply(args[-1], function(f) f[[3]])  # extract RHS of formula
  levels <- levels[!is.na(levels)]
  factor(dplyr::case_when(...), levels=levels)
}

## TTE calculations
tte <- function(origin_date, event_date, censor_date, na.censor=FALSE){
  # returns time-to-event date or time to censor date, which is earlier
  if (na.censor)
    time <- event_date-origin_date
  else
    time <- pmin(event_date-origin_date, censor_date-origin_date, na.rm=TRUE)
  as.numeric(time)
}

## Tidy cox model outputs
tidy_coxph <- function(x, conf.int = TRUE, conf.level = .95, exponentiate = TRUE, ...) {
  
  # to use Wald CIs instead of profile CIs.
  ret <- broom::tidy(x, conf.int = FALSE, conf.level = conf.level, exponentiate = exponentiate)
  
  if(conf.int){
    ci <- confint.default(x, level = conf.level)
    if(exponentiate){ci = exp(ci)}
    ci <- as_tibble(ci, rownames = "term")
    names(ci) <- c("term", "conf.low", "conf.high")
    
    ret <- dplyr::left_join(ret, ci, by = "term")
  }
  ret
}

### Redact vaccine frequency table
redact_table = function(table_input) {
  if (any(table_input$Freq<=100)) {
    # Sum any combinations with less than or equal to 100 records
    other = data.frame(Var1 = "other", Freq = sum(table_input$Freq[table_input$Freq<=100]))
    table_input = rbind(table_input, other)
    table_input = subset(table_input, Freq>100 | Var1=="other")
  }
  
  # Round to nearest 5
  rounded_n = round_any(sum(table_input$Freq),5)
  table_input$Freq = round_any(table_input$Freq,5)
  table_input$Prop = round(table_input$Freq/rounded_n*100,1)
  table_input$Non_Freq = rounded_n - table_input$Freq
  
  # Redact any remaining cell counts <=10
  table_input$Freq[table_input$Freq<=10 | table_input$Non_Freq<=10] = "[Redacted]"
  table_input$Prop[table_input$Freq<=10  | table_input$Non_Freq<=1] = "[Redacted]"
  table_input <- table_input %>% select(-Non_Freq)
  
  names(table_input) = c("Combination", "Frequency", "Percent")
  table_input
}

## Extract coefficients and p values from coxph model output
extract_model = function(model_output) {
  summary = data.frame(
    term = rownames(summary(model_output)[[7]]),
    N = model_output$n, 
    estimate = round(exp(model_output$coefficients),2), 
    conf.low = round(summary(model_output)[[8]][,3],2), 
    conf.high = round(summary(model_output)[[8]][,4],2), 
    p.value = round(summary(model_output)[[7]][,5],4)
  )
}

## Rounding function
round_any = function(x, accuracy, f=round) {
  f(x/accuracy) * accuracy
}
