library(GeRnika)
options(browser = 'zen-browser')

# generate data
inst <- create_instance(
	n = 10,
	m = 4,
	k = 1,
	selection = "neutral")

B <- inst$B

# Create a new 'Phylotree' object
# on the basis of the B matrix
phylotree1 <- B_to_phylotree(B = B)

plot(phylotree1)
