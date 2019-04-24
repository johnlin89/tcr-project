tbl1 = read.table("/Users/linjo/Desktop/tcr-project-desktop/snpFlipResults/windowCovar_flip.1.Q")
barplot(t(as.matrix(tbl1)), col = rainbow(3),
          xlab = "Individual #", ylab = "Ancestry", border = NA)
tbl2 = read.table("/Users/linjo/Desktop/tcr-project-desktop/snpFlipResults/windowCovar_flip.2.Q")
barplot(t(as.matrix(tbl2)), col = rainbow(3),
        xlab = "Individual #", ylab = "Ancestry", border = NA)
tbl3 = read.table("/Users/linjo/Desktop/tcr-project-desktop/snpFlipResults/windowCovar_flip.3.Q")
barplot(t(as.matrix(tbl3)), col = rainbow(3),
        xlab = "Individual #", ylab = "Ancestry", border = NA)
tbl4 = read.table("/Users/linjo/Desktop/tcr-project-desktop/snpFlipResults/windowCovar_flip.4.Q")
barplot(t(as.matrix(tbl4)), col = rainbow(3),
        xlab = "Individual #", ylab = "Ancestry", border = NA)
tbl5 = read.table("/Users/linjo/Desktop/tcr-project-desktop/snpFlipResults/windowCovar_flip.5.Q")
barplot(t(as.matrix(tbl5)), col = rainbow(3),
        xlab = "Individual #", ylab = "Ancestry", border = NA)
