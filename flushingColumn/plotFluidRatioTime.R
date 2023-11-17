# load pyFlushing results and plot surface data in time

# Chiara P. Montagna, INGV Pisa, 4/2022
# chiara.montagna@ingv.it

# ============ TO BE SET ==============================
                                        # select pressures (depths) to plot
pressures<-c(2.e7,1.e8,2.e8)                                   
# =====================================================

setwd('../pyResults/')

                                        # load .out files in a list
columns <- lapply(Sys.glob("*.out"), read.delim, header=FALSE)

                                        # remove headers
colNoHead<-lapply(columns, function(x) x[-1,])                 
colNoHead<-lapply(colNoHead, function(x) x[-1,])

                                        # make data numeric
colNumbers <- lapply(colNoHead, function(x) {
    x <- as.data.frame(sapply(x,as.numeric))
    return(x)
})                                               

                                        # select pressures
selectedRows <- lapply(colNumbers, function(x) {
    x[x$V1 %in% pressures, ]
})

                                        # plot single pressures
pAtm <- lapply(selectedRows, function(x) {
    x[x$V1 %in% 2e7, ]
})

plot(c(0),((1 - pAtm[[1]]$V7)/pAtm[[1]]$V7), xlim = range(c(0,1000)), ylim =
range(c(0,40)), ylab = "Excess CO2/H2O mass", xlab = "Mass of added
fluids")

for (i in 1:length(pAtm)) {
    points(i, ((1 - pAtm[[i]]$V7)/pAtm[[i]]$V7), add=TRUE, cex = 0.01)}

p100 <- lapply(selectedRows, function(x) {
    x[x$V1 %in% 1e8, ]
})

for (i in 1:length(p100)) {
    points(i, ((1 - p100[[i]]$V7)/p100[[i]]$V7), add=TRUE, cex = 0.01, col = "red")}

p200 <- lapply(selectedRows, function(x) {
    x[x$V1 %in% 2e8, ]
})

for (i in 1:length(p200)) {
    points(i, ((1 - p200[[i]]$V7)/p200[[i]]$V7), add=TRUE, cex = 0.01, col = "green")}

legend(0,40, c("Surface","100 MPa","200 MPa"), col=c("black","red","green"), lty=1:2, cex=0.8)


