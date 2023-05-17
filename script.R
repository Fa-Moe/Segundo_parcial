microarreglo <- read.table("de_genes.csv",sep=",", header=T)
microarreglo



####boxplot
control<- microarreglo[,4]
experimental<- microarreglo[,4]*2^microarreglo[,3]
boxplot(control,experimental,ylab="Expresión",xlab="Control                      Experimental")
#Probablemente sea importante aplicar un método de normalización o descartar los outliers en la expresión experimental


###Pre procesamiento de datos
sobre_signif<-ifelse(microarreglo[,3]>0&microarreglo[,6]<0.05,1,0)
sub_signif<-ifelse(microarreglo[,3]<0&microarreglo[,6]<0.05,-1,0)

direccion_signif<- sobre_signif+sub_signif

microarreglo<-cbind(microarreglo,direccion_signif)
head(microarreglo)

#Encuentra los genes que se sobre-expresan

nombres_genes_sobre_expresados <- subset(microarreglo,subset=direccion_signif==1)[,2]
head(nombres_genes_sobre_expresados)
nombres_genes_sub_expresados <- subset(microarreglo,subset=direccion_signif==-1)[,2]
head(nombres_genes_sub_expresados)


write.table(nombres_genes_sobre_expresados,file="nombres_genes_sobre_expresados.csv",sep=",")
####Hacer un volcano plot



any(microarreglo[,6]>0.05) #No hay ningún gen con un P valor mayor a 0.05, todo es significativo

plot(-log10(microarreglo[,6])~microarreglo[,3],col=microarreglo[,8]+3,pch=20,ylab="-log10(p)",xlab="log2FC",ylim=c(0,14))
#Los rojos están sub-expresados significativamente, los azules están sobre-expresados significativamente.