hmdesign2 <-
function (n1, n2) 
{
    n = n1 + n2
    c1 = rep(1, n)
    c2 = 1:n
    c3 = c(rep(0, n1), rep(1, n2))
    c4 = c(rep(0, (n1 + 1)), 1:(n2 - 1))
    hmdesign2 = cbind(c1, c2, c3, c4)
    hmdesign2
}
