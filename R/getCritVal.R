getCritVal <- function(n, alpha = 0.05, filter, r = 1e4, nq = n, options = NULL, stat = NULL, messages = NULL) {
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n)) {
    stop("number of observations 'n' must be a single positive integer")
  }
  
  if (!is.integer(n)) {
    n <- as.integer(n + 1e-6)
  }
  
  if (n < 1L) {
    stop("number of observations 'n' must be a single positive integer")
  }
  
  if (!is.numeric(nq) || length(nq) != 1 || !is.finite(nq)) {
    stop("nq must be a single integer greather or equal than the number of observations 'n'")
  }
  
  if (!is.integer(nq)) {
    nq <- as.integer(nq + 1e-6)
  }
  
  if (nq < n) {
    stop("nq must be a single integer greather or equal than the number of observations 'n'")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a probability, i.e. a single numeric between 0 and 1")
  }
  
  if (!is(filter, "lowpassFilter")) {
    stop("filter must be an object of class 'lowpassFilter'")
  }
  
  stepR::critVal(output = "value", alpha = alpha, n = n, nq = nq, family = "mDependentPS",
                 intervalSystem = "dyaLen", lengths = 2^(0:as.integer(floor(log2(n)))), penalty = "sqrt",
                 stat = stat, filter = filter, r = r, options = options, messages = messages)
}
