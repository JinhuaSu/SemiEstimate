# SemiEstimate

## designer

S3 usage

functional object

there is a global varible s3 -> final result

there are different small function to modify it

function factory + s3

## usage

a function:

build jac_list: check name is correct

```
new_Date <- function(x = double()) {
  stopifnot(is.double(x))
  structure(x, class = "Date")
}

new_Date(c(-1, 0, 1))
#> [1] "1969-12-31" "1970-01-01" "1970-01-02"
```

semislv <- function(theta0, lambda0, Phi_fn, Psi_fn, jac = list(), ...,
method = c("iterative", "implicit"), jacobian=FALSE, control=list())

all build function should be the constructor:

https://adv-r.hadley.nz/s3.html

## S3

---

eqfns(class):

$Phi_fn

$Psi_fn

---

jac(class):

$Phi_der_theta_fn

$Phi_der_lambda_fn

$Psi_der_theta_fn

$Psi_der_lambda_fn

constructor: new_jac -> function()
validator: check the expression name if there is (iter2)

---

quasijac(class):

$Phi_der_theta_fn

$Phi_der_lambda_fn

$Psi_der_theta_fn

$Psi_der_lambda_fn

constructor: new_jac -> function()
validator: check the expression name if there is (iter2)

---

semijac(class):

$Phi_der_theta_fn

$Phi_der_lambda_fn

$Psi_der_theta_fn

$Psi_der_lambda_fn

constructor: new_jac -> function()
validator: check the expression name if there is (iter2)

---

diyjac(class):

$ordered_fn

$itermedials(class)

$return_fn

constructor: new_jac -> function()
validator: check the expression name if there is (iter2)

---

iterspace(class): -> {"ITAT","IPAT","ITHM","IPHM"}

$initials(base list)

$eqfns

$jac_like

$iter_step

$update_delta

$parameters(base list): copy from initial at the step 1

---

resspace(list)

iterspace -> respace

## fn

generic:

update(iterspace) -> (iterspace, iter_over_flag)

update.ITAT

update.IPAT

update.ITHM

update.IPHM

savestats(resspace, iterspace)
