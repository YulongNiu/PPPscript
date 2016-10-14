library('PhyloProfile')
library('ape')

treeFile <- '(((((t1, t2), t3), (t4, t5)), (((t6, t7), (t8, t9)), t10)), ((t11, t12), t13));'
tree <- read.tree(text = treeFile)



treeFile <- '(((t6, t7), (t8, t9)), t10);'
tree <- read.tree(text = treeFile)
(((t1, t2), t3), (t4, t5))
(((t6, t7), (t8, t9)), t10)
((t11, t12), t13)
(((((t1, t2), t3), (t4, t5)), (((t6, t7), (t8, t9)), t10)), ((t11, t12), t13))

tmp1 <- c(1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1)
tmp2 <- c(rep(0:1, each = 5), rep(1, 3))
DolloDist(tree$edge, nodepath(tree), tmp1, tmp1)
