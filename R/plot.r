number_ticks <- function(n) {function(limits) pretty(limits, n)}
plotter <- function(x, index="NRMSE", fileformat="png", xlab) {
    x.transform <- cbind(melt(x, variable.name="Method", value.name=index),xlab)
    if (index == "NRMSE") {
        p1<-ggplot(x.transform,aes(x=xlab,y=NRMSE,colour=Method))+
            geom_line(aes(group=Method),size=1.5)+
            geom_point(aes(shape=Method),size=7)+
            xlab("Missing Rate (%)")+
            ylab("NRMSE")
        if (fileformat == "png") {
            png(filename="nrmse.png",width=1600,height=600,bg="white",res=NA)
            print(p1)
            dev.off()
        } else if (fileformat == "tiff") {
            tiff(filename="nrmse.tiff",width=1600,height=600,bg="white",res=NA)
            print(p1)
            dev.off()
        } else if (fileformat == "jpg") {
            jpeg(filename="nrmse.jpg",width=1600,height=600,bg="white",res=NA)
            print(p1)
            dev.off()
        } else {
            bmp(fileformat == "nrmse.bmp", width=1600, height=600, bg="white", res=NA)
            print(p1)
            dev.off()
        }
    } else if (index == "CPP") {
        p1<-ggplot(x.transform,aes(x=xlab,y=CPP,colour=Method))+
            geom_line(aes(group=Method),size=1.5)+
            geom_point(aes(shape=Method),size=7)+
            scale_x_continuous(breaks=number_ticks(20))+
            xlab("Missing Rate (%)")+
            ylab("CPP")       
        if (fileformat == "png") {
            png(filename="cpp.png",width=1600,height=600,bg="white",res=NA)
            print(p1)
            dev.off()
        } else if (fileformat == "tiff") {
            tiff(filename="cpp.tiff",width=1600,height=600,bg="white",res=NA)
            print(p1)
            dev.off()
        } else if (fileformat == "jpg") {
            jpeg(filename="cpp.jpg",width=1600,height=600,bg="white",res=NA)
            print(p1)
            dev.off()
        } else {
            bmp(fileformat == "cpp.bmp", width=1600, height=600, bg="white", res=NA)
            print(p1)
            dev.off()
        }
    } else {
        p1<-ggplot(x.transform,aes(x=xlab,y=BLCI,colour=Method))+
            geom_line(aes(group=Method),size=1.5)+
            geom_point(aes(shape=Method),size=7)+
            scale_x_continuous(breaks=number_ticks(20))+
            xlab("Missing Rate (%)")+
            ylab("BLCI")
        if (fileformat == "png") {
            png(filename="blci.png",width=1600,height=600,bg="white",res=NA)
            print(p1)
            dev.off()
        } else if (fileformat == "tiff") {
            tiff(filename="blci.tiff",width=1600,height=600,bg="white",res=NA)
            print(p1)
            dev.off()
        } else if (fileformat == "jpg") {
            jpeg(filename="blci.jpg",width=1600,height=600,bg="white",res=NA)
            print(p1)
            dev.off()
        } else {
            bmp(fileformat == "blci.bmp", width=1600, height=600, bg="white", res=NA)
            print(p1)
            dev.off()
        }
    }
}
