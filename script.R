#########Lectura
microarreglo <- read.table("de_genes.csv",sep=",", header=T)
microarreglo



#########boxplot
control<- microarreglo[,4]
experimental<- microarreglo[,4]*2^microarreglo[,3]
boxplot(control,experimental,ylab="Expresión",xlab="Control WT                      Experimental A")
#Probablemente sea importante aplicar un método de normalización o descartar los outliers en la expresión experimental


#########Pre procesamiento de datos
sobre_signif<-ifelse(microarreglo[,3]>0&microarreglo[,6]<0.05,1,0)
sub_signif<-ifelse(microarreglo[,3]<0&microarreglo[,6]<0.05,-1,0)

direccion_signif<- sobre_signif+sub_signif

microarreglo<-cbind(microarreglo,direccion_signif)
head(microarreglo)

########Encuentra los genes que se sobre-expresan y sub-expresan

nombres_genes_sobre_expresados <- subset(microarreglo,subset=direccion_signif==1)[,2]
head(nombres_genes_sobre_expresados)
nombres_genes_sub_expresados <- subset(microarreglo,subset=direccion_signif==-1)[,2]
head(nombres_genes_sub_expresados)


########conteo
num_sobre<-length(nombres_genes_sobre_expresados)
num_sub  <-length(nombres_genes_sub_expresados)
num_sobre
num_sub

num_sobre+num_sub==dim(microarreglo)[1] #Todos los renglones están o sub o sobre-expresados


#######exportar tablas
write.table(nombres_genes_sobre_expresados,file="nombres_genes_sobre_expresados.csv",sep=",",row.names=F)
write.table(nombres_genes_sub_expresados,file="nombres_genes_sub_expresados.csv",sep=",",row.names=F)


######Grafique un PCA
datos_genes<-cbind(control,experimental)#funciona bien

PCA_de_microarreglo <- prcomp(t(datos_genes),center = T,scale. = T)
summary(PCA_de_microarreglo) #el primer componente principal lo puede diferenciar básicamente todo, el segundo no es necesario.
biplot(PCA_de_microarreglo) #Generamos gráfica. La línea roja aparece porque representa los "unscaled axes" de los componentes.
#vemos que el control y el experimental se diferencían bien.

#######Volcano plot

any(microarreglo[,6]>0.05) #No hay ningún gen con un P valor mayor a 0.05, todo es significativo

plot(-log10(microarreglo[,6])~microarreglo[,3],col=microarreglo[,8]+3,pch=20,ylab="-log10(p)",xlab="log2FC",ylim=c(0,14))
#Los rojos están sub-expresados significativamente, los azules están sobre-expresados significativamente. Todos los genes son significativos.



####Los identificadores de los genes no retornaron ningún resultado al analizar con GeneOntology, se adjunta el archivo imagen1, imagen2 para ver 
#que ocurría (no podía encontrar los genes ni en homo sapiens ni en mus musculus) y la imagen3 muestra los identificadores en la página

########Ejercicio 2

BiocManager::install("rtracklayer")
library(Biobase)
library(IRanges)
library(rtracklayer)

# Cuidado! Importar el archivo completo puede requerir ~2-3 GB de RAM.
human = import("Homo_sapiens_exons_small_e70.gtf")

human[,1]
class(human)

seqnames(human)
ranges(human)
strand(human)
mcols(human)
table(mcols(human)$gene_biotype)
mcols(human) = mcols(human)[,c("source","gene_id","gene_name","gene_biotype")]


######Cómo le harían para quedarse exclusivamente con las anotaciones de "miRNA"?
miRNA_humano <- subset(human,subset=mcols(human[,1])[,1]=="miRNA")
miRNA_humano

######y solamente aquellas anotaciones de la cadena "-"?
cadena_neg_humano<-subset(human,subset=strand(human)=="-")
cadena_neg_humano

#La función subset es muy útil, notese que la base utilizada en la condición no tiene que ser igual a la base
#de la que seleccionamos


library(Rsamtools)

what = c("rname", "strand", "pos", "qwidth")
param = ScanBamParam(what=what)

bam = scanBam("human_mapped_small.bam", param=param)

class(bam)

lapply(bam, names)

mapGR = GRanges(
  seqnames = bam[[1]]$rname,
  ranges   = IRanges(start=bam[[1]]$pos, width=bam[[1]]$qwidth),
  strand   = bam[[1]]$strand
)

mapGR

mcols(human)$counts = countOverlaps(human, mapGR)

mcols(human)

typeCounts = aggregate(mcols(human)$counts, by=list("biotype"=mcols(human)$gene_biotype), sum)

typeCounts

geneCounts = aggregate(mcols(human)$counts, by=list("id"=mcols(human)$gene_name), sum)

head(geneCounts)

minCount = 40000
typeCountsHigh = typeCounts[typeCounts$x > minCount,]
typeCountsHigh = typeCountsHigh[order(typeCountsHigh$x),]
typeCountsHigh = rbind(data.frame("biotype"="other",
                                  "x"=sum(typeCounts$x[typeCounts$x <= minCount])),
                       typeCountsHigh)

pie(typeCountsHigh$x, labels=typeCountsHigh$biotype, col=rev(rainbow(nrow(typeCountsHigh))),
    main="Number of aligned reads per biotype")
####Aquí podemos ver que hay más lecturas asociadas a miRNA que a protein-coding genes y una porción insignificante
#o nula de otros