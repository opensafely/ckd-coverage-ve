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
    time <- date(event_date)-date(origin_date)
  else
    time <- pmin(date(event_date)-date(origin_date), date(censor_date)-date(origin_date), na.rm=TRUE)
  as.numeric(time)
}

censor_indicator <- function(event_date, censor_date){
  # returns 0 if event_date is censored by censor_date, or if event_date is NA. Otherwise 1
  dplyr::if_else((event_date>censor_date) | is.na(event_date), FALSE, TRUE)
}

rrCI_exact <- function(n, pt, ref_n, ref_pt, accuracy=0.001){
  
  # use exact methods if incidence is very low for immediate post-vaccine outcomes
  
  rate <- n/pt
  ref_rate <- ref_n/ref_pt
  rr <- rate/ref_rate
  
  ll = ref_pt/pt * (n/(ref_n+1)) * 1/qf(2*(ref_n+1), 2*n, p = 0.05/2, lower.tail = FALSE)
  ul = ref_pt/pt * ((n+1)/ref_n) * qf(2*(n+1), 2*ref_n, p = 0.05/2, lower.tail = FALSE)
  
  paste0("(", scales::number_format(accuracy=accuracy)(ll), "-", scales::number_format(accuracy=accuracy)(ul), ")")
  
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
  table_input$Prop[table_input$Freq<=10  | table_input$Non_Freq<=10] = "[Redacted]"
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

## Extract coefficients and p values from logistic model output (quick method)
extract_model_logistic = function(variable_name, model_output) {
  summary = data.frame(
    term = paste0(as.character(unlist(model_output$xlevels))),
    N = length(model_output$fitted.values), 
    estimate = round(exp(model_output$coefficients),2), 
    conf.low = round(exp(confint.default(model_output))[,1],2), 
    conf.high = round(exp(confint.default(model_output))[,2],2), 
    p.value = round(summary(model_output)$coefficients[,4],5)
  )
  # Set intercept values to NA
  summary[1,3:6] = NA
  return(summary)
}

## Rounding function
round_any = function(x, accuracy, f=round) {
  f(x/accuracy) * accuracy
}


timesince_cut <- function(time_since, breaks, prelabel="pre", prefix=""){
  
  # this function defines post-vaccination time-periods at `time_since`,
  # delimited by `breaks`
  
  # note, intervals are open on the left and closed on the right
  # so at the exact time point the vaccination occurred, it will be classed as "pre-dose".
  
  time_since <- as.numeric(time_since)
  time_since <- if_else(!is.na(time_since), time_since, Inf)
  
  breaks_aug <- unique(c(-Inf, breaks, Inf))
  
  lab_left <- breaks+1 
  lab_right <- lead(breaks) 
  label <- paste0(lab_left, "-", lab_right)
  label <- str_replace(label,"-NA", "+")
  labels <- paste0(prefix, c(prelabel, label))
  
  period <- cut(time_since, breaks=breaks_aug, labels=labels, include.lowest=TRUE)
  period
}

timesince_cut_end <- function(time_since, breaks, prefix=""){

  # this function defines post-vaccination time-periods at `time_since`,
  # delimited by `breaks`

  # note, intervals are open on the left and closed on the right
  # so at the exact time point the vaccination occurred, it will be classed as "pre-dose".

  stopifnot("time_since should be strictly non-negative" = time_since>=0)
  time_since <- as.numeric(time_since)
  time_since <- if_else(!is.na(time_since), time_since, Inf)

  lab_left <- breaks[-length(breaks)]+1
  lab_right <- breaks[-1]
  label <- paste0(lab_left, "-", lab_right)
  labels <- paste0(prefix, label)

  #labels0 <- cut(c(breaks, Inf), breaks_aug)
  #labels <- paste0(prefix, c(prelabel, as.character(labels0[-1])))
  period <- cut(time_since, breaks=breaks, labels=labels, include.lowest=FALSE)

  period
}

redactor2 <- function(n, threshold=5, x=NULL){
  
  # given a vector of frequencies (n), this returns a redacted vector (if x is NULL) or
  # reaction of a secondary vector based on frequencies in the first (if x is not nULL).
  # using the following rules:
  # a) the frequency is <= the redaction threshold and
  # b) if the sum of redacted frequencies in a) is still <= the threshold, then the
  # next largest frequency is also redacted
  
  
  stopifnot("n must be non-missing" = any(!is.na(n)))
  stopifnot("n must non-negative" = any(n>=0))
  stopifnot("n must non-negative" = any(n>=0))
  
  if(is.null(x)){
    x <- n
  }
  
  if(!is.null(x)){
    stopifnot("x must be same length as n" = length(n) == length(x))
  }
  
  
  
  n <- as.integer(n)
  leq_threshold <- dplyr::between(n, 1, threshold)
  n_sum <- sum(n)
  
  # redact if n is less than or equal to redaction threshold
  redact <- leq_threshold
  
  # also redact next smallest n if sum of redacted n is still less than or equal to threshold
  if((sum(n*leq_threshold) <= threshold) & any(leq_threshold)){
    redact[which.min(dplyr::if_else(leq_threshold, n_sum+1L, n))] = TRUE
  }
  
  
  typedNA <- NA
  mode(typedNA) <- typeof(x)
  
  redacted <- dplyr::if_else(redact, typedNA, x)
  
  redacted
}


## PLR function from HCW comparative effectiveness analysis
sample_nonoutcomes_n <- function(had_outcome, id, n){
  # TRUE if outcome occurs,
  # TRUE with probability of `prop` if outcome does not occur
  # FALSE with probability `prop` if outcome does occur
  # based on `id` to ensure consistency of samples
  
  # `had_outcome` is a boolean indicating if the subject has experienced the outcome or not
  # `id` is a identifier with the following properties:
  # - a) consistent between cohort extracts
  # - b) unique
  # - c) completely randomly assigned (no correlation with practice ID, age, registration date, etc etc) which should be true as based on hash of true IDs
  # - d) is an integer strictly greater than zero
  # `proportion` is the proportion of nonoutcome patients to be sampled
  
  (dplyr::dense_rank(dplyr::if_else(had_outcome, 0L, id)) - 1L) <= n
  
}

## PLR function from HCW comparative effectiveness analysis
sample_weights <- function(had_outcome, sampled){
  # `had_outcome` is a boolean indicating if the subject has experienced the outcome or not
  # `sampled` is a boolean indicating if the patient is to be sampled or not
  case_when(
    had_outcome ~ 1,
    !had_outcome & !sampled ~ 0,
    !had_outcome & sampled ~ sum(!had_outcome)/sum((sampled) & !had_outcome),
    TRUE ~ NA_real_
  )
}

glance_plr <- function(model){
  tibble(
    AIC=model$aic,
    df.null=model$df.null,
    df.residual=model$df.residual,
    deviance=model$deviance,
    null.deviance=model$null.deviance,
    nobs=length(model$y)
  )
}

tidy_plr <- function(model, conf.int=TRUE, conf.level=0.95, exponentiate=FALSE, cluster){
  
  # create tidy dataframe for coefficients of pooled logistic regression
  # using robust standard errors
  robustSEs <- lmtest::coeftest(model, vcov. = sandwich::vcovCL(model, cluster = cluster, type = "HC0")) %>% broom::tidy(conf.int=FALSE, exponentiate=exponentiate)
  robustCIs <- lmtest::coefci(model, vcov. = sandwich::vcovCL(model, cluster = cluster, type = "HC0"), level = conf.level) %>% tibble::as_tibble(rownames="term")
  robust <- dplyr::inner_join(robustSEs, robustCIs, by="term")
  
  robust %>%
    rename(
      conf.low=`2.5 %`,
      conf.high=`97.5 %`
    ) %>%
    mutate(
      or = exp(estimate),
      or.ll = exp(conf.low),
      or.ul = exp(conf.high),
    )
  
}
