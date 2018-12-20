# Acknowledgement: This file is modified from gWQS.


# functtion to define the objective function
objfn <- function(initp, Q, y, family, covrts, b1_constr){

  b0 = initp[1]
  b1 = initp[2]
  w = initp[3:(dim(Q)[2]+2)]
  cvrt_coeff = initp[(dim(Q)[2]+3):length(initp)]

  if(dim(covrts)[2] == 0) term = b0 + b1*Q%*%w
  else term = b0 + b1*Q%*%w + covrts%*%cvrt_coeff

  if (family == "gaussian") f = sum((y - term)^2)
  else if (family == "binomial") {
    p = 1/(1 + exp(-term))
    f = -sum(y * log(p) + (1 - y) * log(1 - p))
  }

  return(f)
}


# function to define the equality constraint
linconst = function(initp, Q, y, family, covrts, b1_constr){

  wsum=sum(initp[3:(dim(Q)[2] + 2)])

  return(wsum)
}


# function to determine bounded prameters
bounded_param = function(initp, Q, y, family, covrts, b1_constr){

  if (b1_constr) b_p = initp[2:(dim(Q)[2] + 2)]
  else b_p = initp[3:(dim(Q)[2] + 2)]

  return(b_p)
}


# function to define the lower bounds
LBound = function(Q, b1_pos, b1_constr){

  if (b1_constr) {
    b1l = ifelse(b1_pos, 0, -Inf)
    LB = c(b1l, rep(0, dim(Q)[2]))
  }
  else LB = rep(0, dim(Q)[2])

  return(LB)
}


# function to define the upper bounds
UBound = function(Q, b1_pos, b1_constr){

  if (b1_constr) {
    b1u = ifelse(b1_pos, Inf, 0)
    UB = c(b1u, rep(1, dim(Q)[2]))
  }
  else UB = rep(1, dim(Q)[2])

  return(UB)
}


# function to define the parameters initial values
values.0 = function(lenp, y, covrts, b1_pos, family){

  y = as.matrix(y)
  b1_0 = ifelse(b1_pos, 0.001, -0.001)

  if (dim(covrts)[2] == 0){
    val.0 = c(0, b1_0, rep(1/lenp, lenp))
    names(val.0) = c("b0", "b1", paste0("w", 1:lenp))
  }
  else {
    if (family == "gaussian") fit = lm(y ~ covrts)
    else if (family == "binomial") fit = glm(y ~ covrts, family = binomial(link = "logit"))
    bj = coef(fit)[-1]
    val.0 = c(0, b1_0, rep(1/lenp, lenp), bj)
    names(val.0) = c("b0", "b1", paste0("w", 1:lenp), paste0("b", 2:(dim(covrts)[2]+1)))
  }

  return(val.0)
}


# optimization function to estimate the weights
optim.f <- function(Q, y, b1_pos, b1_constr, family, covrts){
  Q = as.matrix(Q)
  if(dim(covrts)[2] > 0) covrts = as.matrix(covrts)
  lenp = dim(Q)[2]
  initp = values.0(lenp, y, covrts, b1_pos, family)
  LowB = LBound(Q, b1_pos, b1_constr)
  UpB = UBound(Q, b1_pos, b1_constr)

  opt_res = tryCatch(Rsolnp::solnp(pars = initp, fun = objfn, eqfun = linconst, eqB = 1,
                           ineqfun = bounded_param, ineqLB = LowB, ineqUB = UpB,
                           control = list(trace = 0), Q = Q, y = y, family = family,
                           covrts = covrts, b1_constr = b1_constr),
                     error = function(e) NULL)
  if(!is.null(opt_res)) {
    par_opt = opt_res$pars
    conv = opt_res$convergence
    nfuneval = opt_res$nfuneval
  } else {par_opt=initp; conv=2; nfuneval=0}
  out = list(par_opt, conv, nfuneval)
  names(out) = c("par_opt", "conv", "nfuneval")

  return(out)
}


# function that fit the wqs model
model.fit <- function(Q, y, wghts, family, covrts, wqs2){

  Q = as.matrix(Q)
  wqs = Q%*%wghts
  wqs = as.data.frame(wqs)
  names(wqs) = "wqs"

  new_data = cbind(as.data.frame(y), wqs)
  names(new_data)[1] = "y"

  if (dim(covrts)[2] > 0) {
    new_data = cbind(new_data, covrts)
    names(new_data)[c((dim(wqs)[2] + 2):dim(new_data)[2])] = names(covrts)
  }

  if (family == "gaussian") m_f = lm(y ~ ., data = new_data)
  else if (family == "binomial") m_f = glm(y ~ ., data = new_data, family = binomial(link = "logit"))

  mf_out = list(wqs, m_f)
  names(mf_out) = c("wqs", "m_f")

  if (wqs2 == TRUE) {
    wqs_2 = wqs^2
    wqs_2 = as.data.frame(wqs_2)
    names(wqs_2) = "wqs_2"
    wqs = cbind(wqs, wqs_2)
    new_data = cbind(new_data, wqs_2)

    if (family == "gaussian") m_f2 = lm(y ~ ., data = new_data)
    else if (family == "binomial") m_f2 = glm(y ~ ., data = new_data, family = binomial(link = "logit"))

    aov = anova(m_f, m_f2, test = "Chisq")

    mf_out = list(wqs[, 1], m_f, m_f2, aov)
    names(mf_out) = c("wqs", "m_f", "m_f2", "aov")
  }

  return(mf_out)
}


# function that calls the optimization function and the function to fit the model for each bootstrap sample
par.modl.est <- function(data, y_name, q_name, cov_name, b, b1_pos, b1_constr, family, seed){

  index_b = vector("list", b)
  param = vector("list", b)
  b1 = rep(NA, b)
  wght_matrix = matrix(NA, b, length(q_name))
  colnames(wght_matrix) = paste("w", 1:length(q_name), sep = "_")
  conv = rep(NA, b)
  p_val = rep(NA, b)
  nfuneval = rep(NA, b)

  if (family == "gaussian") ts = "t"
  else if (family == "binomial") ts = "z"

  for (i in 1:b) {

    # create bootstrap sample
    index_b[[i]] = sample(1:nrow(data), nrow(data), replace=TRUE)
    data_b <- data[index_b[[i]],]

    # calculate parameters
    param[[i]] = optim.f(data_b[, q_name, drop = FALSE], data_b[, y_name, drop = FALSE],
                         b1_pos, b1_constr, family, data_b[, cov_name, drop = FALSE])

    # fit the wqs model for each bootstrap sample
    b_fit = model.fit(data_b[, q_name, drop = FALSE], data_b[, y_name, drop = FALSE],
                      param[[i]]$par_opt[3:(2 + length(q_name))], family,
                      data_b[, cov_name, drop = FALSE], wqs2 = FALSE)

    wght_matrix[i,] = param[[i]]$par_opt[3:(2 + length(q_name))]
    b1[i] = b_fit$m_f$coefficients[row.names = "wqs"]
    conv[i] = param[[i]]$conv
    p_val[i] = summary(b_fit$m_f)$coefficients["wqs", gsub("x", ts, "Pr(>|x|)")]
    nfuneval[i] = param[[i]]$nfuneval
  }

  n_non_conv = sum(conv == 2)
  print(paste0("The optimization function did not converge ", n_non_conv, " times"))

  par_model_out = list(wght_matrix, b1, conv, p_val, index_b, nfuneval)
  names(par_model_out) = c("wght_matrix", "b1", "conv", "p_val", "index_b", "nfuneval")

  return(par_model_out)
}
