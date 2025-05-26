#cat BLAST_knob180_W23-DfK.out BLAST_knob180_W23-DfK_R2.out >> BLAST_knob180_W23-DfK_ALL.out

#cat BLAST_knob180_W23-DfL.out BLAST_knob180_W23-DfL_R2.out >> BLAST_knob180_W23-DfL_ALL.out

#cat BLAST_knob180_W23-Ab10.out BLAST_knob180_W23-Ab10_R2.out >> BLAST_knob180_W23-Ab10_ALL.out

setwd("")

AB10 <- read.delim("BLAST_knob180_W23-Ab10_ALL.merge.bed")

K <- read.delim("BLAST_knob180_W23-DfK_ALL.merge.bed")

L <- read.delim("BLAST_knob180_W23-DfL_ALL.merge.bed")

names(AB10) <- c("qseqid", "qstart", "qend")

names(L) <- c("qseqid", "qstart", "qend")

names(K) <- c("qseqid", "qstart", "qend")

AB10$length <- abs(AB10$qend-AB10$qstart)
L$length <- abs(L$qend-L$qstart)
K$length <- abs(K$qend-K$qstart)

AB10_filt <- subset(AB10, length >=30)

#y values are average coverage of the genome determined by samtools

#Ab10 
x=sum(AB10_filt$length)
y=4.60642
#normalized by coverage
AB10_MB=(x/y)/1000000
print(AB10_MB)
#AB10 knob is 91.67075

L_filt <- subset(L, length >=30)

#L
a=sum(L_filt$length)
b=5.89151
L_MB=(a/b)/1000000
print(L_MB)
#DfL knob is 91.41426

K_filt <- subset(K, length >=30)

#K
c=sum(K_filt$length)
d=5.20558
K_MB=(c/d)/1000000
print(K_MB)
#DfK knob is 61.0005

AB10_KNOB = AB10_MB - K_MB
print(AB10_KNOB)
#The knob on Ab10 is 30.67025MB

DFL_KNOB = AB10_MB - L_MB
print(DFL_KNOB)
#The knob on Df-L is only 0.2564874 smaller than the knob on Ab10 
