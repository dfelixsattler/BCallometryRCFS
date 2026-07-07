# Height-diameter prediction equations: Naslund (1937), Curtis (1967), logistic,
# Korf (1939), Weibull, and Richards (1959), plus the unified dispatcher
# ht_from_dbh().  All functions are exported.

#' Predict tree height from DBH using the Naslund (1937) equation
#'
#' @description
#' Computes predicted total tree height from diameter at breast height (DBH)
#' using the Naslund (1937) height-diameter equation:
#' \deqn{H = 1.3 + \frac{DBH^2}{(a + b \cdot DBH)^2}}
#'
#' @param dbh numeric vector. Diameter at breast height (cm).
#' @param a   numeric. Fitted coefficient \eqn{a > 0}.
#' @param b   numeric. Fitted coefficient \eqn{b > 0}.
#'
#' @return numeric vector of predicted heights (m), same length as \code{dbh}.
#' @export
#' @references
#'   Naslund, M. (1937). Skogsforsokanstaltens gallringsförsök i tallskog.
#'   \emph{Meddelanden fran Statens Skogsförsöksanstalt}, 29, 1–169.
#' @seealso \code{\link{ht_from_dbh}}
#' @examples
#' ht_naslund(dbh = c(10, 20, 30), a = 1.5, b = 0.05)
ht_naslund <- function(dbh, a, b) {
  1.3 + (dbh^2) / (a + b * dbh)^2
}


#' Predict tree height from DBH using the Curtis (1967) equation
#'
#' @description
#' Computes predicted total tree height from DBH using the Curtis (1967)
#' height-diameter equation:
#' \deqn{H = 1.3 + a \cdot \left(\frac{DBH}{1 + DBH}\right)^b}
#'
#' @param dbh numeric vector. Diameter at breast height (cm).
#' @param a   numeric. Fitted coefficient \eqn{a > 0} (controls asymptote).
#' @param b   numeric. Fitted coefficient \eqn{b > 0} (controls shape).
#'
#' @return numeric vector of predicted heights (m), same length as \code{dbh}.
#' @export
#' @references
#'   Curtis, R.O. (1967). Height-diameter and height-diameter-age equations for
#'   second-growth Douglas-fir. \emph{Forest Science}, 13(4), 365–375.
#' @seealso \code{\link{ht_from_dbh}}
#' @examples
#' ht_curtis(dbh = c(10, 20, 30), a = 35, b = 0.9)
ht_curtis <- function(dbh, a, b) {
  1.3 + a * (dbh / (1 + dbh))^b
}


#' Predict tree height from DBH using the logistic equation
#'
#' @description
#' Computes predicted total tree height from DBH using a three-parameter
#' logistic height-diameter equation:
#' \deqn{H = 1.3 + \frac{a}{1 + b \cdot e^{-c \cdot DBH}}}
#'
#' @param dbh numeric vector. Diameter at breast height (cm).
#' @param a   numeric. Fitted coefficient \eqn{a > 0} (upper asymptote above 1.3 m).
#' @param b   numeric. Fitted coefficient \eqn{b > 0}.
#' @param c   numeric. Fitted coefficient \eqn{c > 0} (growth rate).
#'
#' @return numeric vector of predicted heights (m), same length as \code{dbh}.
#' @export
#' @seealso \code{\link{ht_from_dbh}}
#' @examples
#' ht_logistic(dbh = c(10, 20, 30), a = 40, b = 5, c = 0.1)
ht_logistic <- function(dbh, a, b, c) {
  1.3 + a / (1 + b * exp(-c * dbh))
}


#' Predict tree height from DBH using the Korf equation
#'
#' @description
#' Computes predicted total tree height from DBH using the Korf (1939)
#' height-diameter equation:
#' \deqn{H = 1.3 + a \cdot e^{-b \cdot DBH^{-c}}}
#'
#' @param dbh numeric vector. Diameter at breast height (cm).
#' @param a   numeric. Fitted coefficient \eqn{a > 0} (upper asymptote above 1.3 m).
#' @param b   numeric. Fitted coefficient \eqn{b > 0}.
#' @param c   numeric. Fitted coefficient \eqn{c > 0}.
#'
#' @return numeric vector of predicted heights (m), same length as \code{dbh}.
#' @export
#' @references
#'   Korf, V. (1939). Príspevek k matematické definici vzrůstového zákona
#'   lesních porostu. \emph{Lesnická Práce}, 18, 339–356.
#' @seealso \code{\link{ht_from_dbh}}
#' @examples
#' ht_korf(dbh = c(10, 20, 30), a = 40, b = 3, c = 0.5)
ht_korf <- function(dbh, a, b, c) {
  1.3 + a * exp(-b * dbh^(-c))
}


#' Predict tree height from DBH using the Weibull equation
#'
#' @description
#' Computes predicted total tree height from DBH using a three-parameter
#' Weibull height-diameter equation:
#' \deqn{H = 1.3 + a \cdot (1 - e^{-b \cdot DBH^c})}
#'
#' @param dbh numeric vector. Diameter at breast height (cm).
#' @param a   numeric. Fitted coefficient \eqn{a > 0} (upper asymptote above 1.3 m).
#' @param b   numeric. Fitted coefficient \eqn{b > 0}.
#' @param c   numeric. Fitted coefficient \eqn{c > 0}.
#'
#' @return numeric vector of predicted heights (m), same length as \code{dbh}.
#' @export
#' @seealso \code{\link{ht_from_dbh}}
#' @examples
#' ht_weibull(dbh = c(10, 20, 30), a = 40, b = 0.05, c = 1.2)
ht_weibull <- function(dbh, a, b, c) {
  1.3 + a * (1 - exp(-b * dbh^c))
}


#' Predict tree height from DBH using the Richards equation
#'
#' @description
#' Computes predicted total tree height from DBH using the Richards (1959)
#' height-diameter equation:
#' \deqn{H = 1.3 + a \cdot (1 - e^{-b \cdot DBH})^c}
#'
#' @param dbh numeric vector. Diameter at breast height (cm).
#' @param a   numeric. Fitted coefficient \eqn{a > 0} (upper asymptote above 1.3 m).
#' @param b   numeric. Fitted coefficient \eqn{b > 0}.
#' @param c   numeric. Fitted coefficient \eqn{c > 0}.
#'
#' @return numeric vector of predicted heights (m), same length as \code{dbh}.
#' @export
#' @references
#'   Richards, F.J. (1959). A flexible growth function for empirical use.
#'   \emph{Journal of Experimental Botany}, 10(29), 290–300.
#' @seealso \code{\link{ht_from_dbh}}
#' @examples
#' ht_richards(dbh = c(10, 20, 30), a = 40, b = 0.05, c = 2)
ht_richards <- function(dbh, a, b, c) {
  1.3 + a * (1 - exp(-b * dbh))^c
}


#' Predict tree height from DBH using a named model form
#'
#' @description
#' A unified dispatcher for all supported height-diameter equation forms.
#' Supply the model name and fitted coefficients; the appropriate equation
#' function is called internally.
#'
#' Supported model forms (all share the 1.3 m breast-height convention):
#' \describe{
#'   \item{naslund}{Naslund (1937): \eqn{H = 1.3 + DBH^2 / (a + b \cdot DBH)^2}}
#'   \item{curtis}{Curtis (1967): \eqn{H = 1.3 + a \cdot (DBH / (1 + DBH))^b}}
#'   \item{logistic}{\eqn{H = 1.3 + a / (1 + b \cdot e^{-c \cdot DBH})}}
#'   \item{korf}{Korf (1939): \eqn{H = 1.3 + a \cdot e^{-b \cdot DBH^{-c}}}}
#'   \item{weibull}{\eqn{H = 1.3 + a \cdot (1 - e^{-b \cdot DBH^c})}}
#'   \item{richards}{Richards (1959): \eqn{H = 1.3 + a \cdot (1 - e^{-b \cdot DBH})^c}}
#' }
#'
#' @param dbh   numeric vector. Diameter at breast height (cm).
#' @param model character. One of \code{"naslund"}, \code{"curtis"},
#'   \code{"logistic"}, \code{"korf"}, \code{"weibull"}, \code{"richards"}.
#'   Case-insensitive.
#' @param a     numeric. Fitted coefficient \eqn{a}.
#' @param b     numeric. Fitted coefficient \eqn{b}.
#' @param c     numeric or \code{NULL}. Fitted coefficient \eqn{c}. Required for
#'   three-parameter models (\code{logistic}, \code{korf}, \code{weibull},
#'   \code{richards}); ignored for two-parameter models (\code{naslund},
#'   \code{curtis}).
#'
#' @return numeric vector of predicted heights (m), same length as \code{dbh}.
#' @export
#' @seealso \code{\link{ht_naslund}}, \code{\link{ht_curtis}},
#'   \code{\link{ht_logistic}}, \code{\link{ht_korf}},
#'   \code{\link{ht_weibull}}, \code{\link{ht_richards}},
#'   \code{\link{fit_hd_model}}
#' @examples
#' # Two-parameter model
#' ht_from_dbh(dbh = c(10, 20, 30), model = "naslund", a = 1.5, b = 0.05)
#'
#' # Three-parameter model
#' ht_from_dbh(dbh = c(10, 20, 30), model = "korf", a = 40, b = 3, c = 0.5)
ht_from_dbh <- function(dbh, model, a, b, c = NULL) {
  model <- tolower(model)
  switch(model,
    naslund  = ht_naslund(dbh, a, b),
    curtis   = ht_curtis(dbh, a, b),
    logistic = ht_logistic(dbh, a, b, c),
    korf     = ht_korf(dbh, a, b, c),
    weibull  = ht_weibull(dbh, a, b, c),
    richards = ht_richards(dbh, a, b, c),
    stop(
      "Unknown model '", model, "'. ",
      "Must be one of: naslund, curtis, logistic, korf, weibull, richards."
    )
  )
}
