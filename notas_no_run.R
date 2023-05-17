##Volcano

plot(exp_tmem_proc[,12]~exp_tmem_proc[,10],col=exp_tmem_proc[,13]+3,pch=20,ylab="-log10(p, T-test)",xlab="log2FC")

##GEO2R

##consultar tarea 3


###Simulación
#################
---
  title: "Medicion expresion a mano"
author: "Farell & Oscar"
date: "`r Sys.Date()`"
output: html_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Primero cargamos la libreria de un fasta de referencia "Acinetobacter baumannii strain 4300STDY7045842 genome assembly, contig: ERS1930315SCcontig000005, whole genome shotgun sequence" (El fasta se descargó de NCBI)</p>
  
  ```{r import-data1, cache=TRUE, fig.width = 10, fig.asp = 1}
library(stringi)

library(seqinr)

secuencia <- read.fasta("/Documentos/Escuela/VI_Genomica_Roberto/Local_genomica/clase/24-abr/sequence.fasta", as.string=T,set.attributes = F,strip.desc = T)


```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Trabajaremos con las primeras 1,500 bases de la secuencia</p>
  
  ```{r import-data2, cache=TRUE, fig.width = 10, fig.asp = 1}

param_tamaño_secuencia<-1500
param_espacio_inicial_min <-0
param_espacio_inicial_max <-25
param_min_tamaño_gen_anotado <-100
param_max_tamaño_gen_anotado <-200
param_numero_lecturas <-3500
param_min_tamaño_lectura <-10
param_max_tamaño_lectura <- 15

Secuencia_Ensamblada <- stri_sub(unlist(secuencia),0,param_tamaño_secuencia)
Secuencia_Ensamblada


```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Para simular una anotación sobre nuestra secuencia, generaremos genes de entre 100 a 200 bases. Para que no se sobrelape, el inicio de cada uno es el final del último gen generado. Al final se generan las posiciones de nuestros "genes" anotados.</p>
  
  ```{r import-data3, cache=TRUE, fig.width = 10, fig.asp = 1}

anotacion_pos<-c(0,0)
final<-0
inicio <- floor(runif(1, min=param_espacio_inicial_min,max=param_espacio_inicial_max))
while (final <= param_tamaño_secuencia){
  
  final <- floor(inicio+runif(1, min=param_min_tamaño_gen_anotado,max=param_max_tamaño_gen_anotado))
  anotacion_pos<-rbind(anotacion_pos,c(inicio,final))
  inicio <- final
}
anotacion_pos<- anotacion_pos[-c(1,dim(anotacion_pos)[1]),]
anotacion_pos

```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Obtenemos con stri_sub ahora los nucleótidos correspondientes a los genes.</p>
  
  ```{r import-data4, cache=TRUE, fig.width = 10, fig.asp = 1}

anotacion_nuc <- NULL
for (i in 1:dim(anotacion_pos)[1]){
  anotacion_nuc <- list(unlist(anotacion_nuc),stri_sub(Secuencia_Ensamblada,anotacion_pos[i,1],anotacion_pos[i,2]))
}
anotacion_nuc <- list(unlist(anotacion_nuc))
anotacion_nuc

```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Ahora generamos un montón de lecturas. Su inicio y longitud son aleatorios.</p>
  
  ```{r import-data5, cache=TRUE, fig.width = 10, fig.asp = 1}

lecturas <- c(0,0)


for(i in 1:param_numero_lecturas){
  
  start <- ceiling(runif(1, min = anotacion_pos[1,1], max = anotacion_pos[dim(anotacion_pos)[1],2]-param_max_tamaño_lectura))
  longitude <- ceiling(runif(1, min = param_min_tamaño_lectura, max = param_max_tamaño_lectura))
  
  lecturas <- rbind(lecturas,c(start,longitude))
}
lecturas <-lecturas[2:dim(lecturas)[1],]

head(lecturas)

```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Trasladamos las lecturas de números a los conjuntos de bases que representa cada lectura.</p>
  
  ```{r import-data6, cache=TRUE, fig.width = 10, fig.asp = 1}

Lecturas_nucleotido <- NULL
for (i in 1:dim(lecturas)[1]){
  Lecturas_nucleotido <- list(unlist(Lecturas_nucleotido),stri_sub(Secuencia_Ensamblada,lecturas[i,1],lecturas[i,1]+lecturas[i,2]))
}
Lecturas_nucleotido <- list(unlist(Lecturas_nucleotido))
head(Lecturas_nucleotido[[1]])

```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Primero convertimos la secuencia base a un formato Biostrings. Luego especificamos los "genes" o la parte anotada con los índices correspondientes, generando un objeto "XStringViews". Al seleccionar con corchetes una de estas views, se genera como si la secuencia fuera un objeto "XString", que se puede usar luego para contar cuantas veces aparece dada una secuencia</p>
  
  ```{r import-data7, cache=TRUE, fig.width = 10, fig.asp = 1}


library(Biostrings)

genoma_2<-DNAString(Secuencia_Ensamblada)
genoma_2

genes_anotados_2<-Views(genoma_2,start=anotacion_pos[,1] ,end=anotacion_pos[,2])
genes_anotados_2[[1]]

```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Luego generamos una matriz en blanco, con tantas columnas como "genes" y tantos renglones como "lecturas". En un bucle, contamos cuantas veces aparece cada lectura en todos los genes. El resultado es una matriz que indica que lectura aparece en que gen</p>
  
  ```{r import-data8, cache=TRUE, fig.width = 10, fig.asp = 1}


matriz_apariciones <- matrix(0,nrow=length(Lecturas_nucleotido[[1]]),ncol=length(genes_anotados_2))

for (i in 1:length(Lecturas_nucleotido[[1]])){
  
  for (j in 1:length(genes_anotados_2)){
    matriz_apariciones[i,j]<-countPattern(Lecturas_nucleotido[[1]][i],genes_anotados_2[[j]])
  }
}
head(matriz_apariciones)


```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Para resumir la información, generamos un renglón que nos diga cuantas lecturas totales existen por gen.</p>
  
  ```{r import-data9, cache=TRUE, fig.width = 10, fig.asp = 1}
resumen_apariciones <- matrix(0,ncol=length(genes_anotados_2))

for (i in 1:dim(resumen_apariciones)[2]){
  resumen_apariciones[1,i]<-sum(matriz_apariciones[,i])
}

resumen_apariciones

```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Y ya podemos graficar cuantas lecturas hubo en cada gen.</p>
  
  ```{r import-data10, cache=TRUE, fig.width = 10, fig.asp = 1}

plot(resumen_apariciones[1,],xlab="indice(#) del gen",ylab="lecturas")

```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Luego rescatamos el tamaño de cada gen. Se podría sacar a partir de anotación_pos con una operación pero este método que utiliza el tamaño del gen ya leído como ADN también es útil.</p>
  
  ```{r import-data11, cache=TRUE, fig.width = 10, fig.asp = 1}

tamaño_gen_store <- c()
tamaño_gen_act <- NULL
for (i in 1:length(genes_anotados_2)){
  tamaño_gen_act <- length(genes_anotados_2[[i]])
  tamaño_gen_store <- cbind(tamaño_gen_store,tamaño_gen_act)
}
head(tamaño_gen_store) 

```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Podemos graficar el tamaño del gen respectivo a su longitud, para ver como se relacionan estos. Como la gráfica parece indicar una relación lineal, sacamos el índice de correlación y efectivamente vemos que hay una muy fuerte correlación entre el tamaño del gen y el número de lecturas que genera.</p>
  
  ```{r import-data12, cache=TRUE, fig.width = 10, fig.asp = 1}
plot(resumen_apariciones[1,]~tamaño_gen_store[1,],ylab="Número de lecturas en ese gen",xlab="Longitud del gen")

cor(resumen_apariciones[1,],tamaño_gen_store[1,])

```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  También está la opción de utilizar GGally para visualizar la correlación entre las variables.</p>
  
  ```{r import-data13, cache=TRUE, fig.width = 10, fig.asp = 1}

library(ggplot2)
library(GGally)
lecturas_tamaño<-as.data.frame(unname(t(rbind(resumen_apariciones[1,],tamaño_gen_store[1,]))))
colnames(lecturas_tamaño) <- c("lecturas","tamaño")
ggpairs(lecturas_tamaño)

```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Luego podemos normalizar las lecturas por gen para tomar en cuenta su tamaño. Vemos que los valores de expresión entre los genes son parecidos a excepción de un poco de ruido. Esto fue debido a que se utilizó runif (aleatorio uniforme) para generar las lecturas, por lo que cada gen tendría el mismo número de transcritos una vez se toma en cuenta su tamaño. Esperaríamos que conforme las lecturas crezcan, la diferencia relativa entre los valores de expresión disminuya.</p>
  
  ```{r import-data14, cache=TRUE, fig.width = 10, fig.asp = 1}

expresion<- unname(resumen_apariciones[1,]/tamaño_gen_store[1,])
expresion

```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Generamos una tabla de abundancias</p>
  
  
  ```{r import-data14b, cache=TRUE, fig.width = 10, fig.asp = 1}

tabla_de_abundancia<- t(unname(rbind(as.numeric(resumen_apariciones[1,]),paste0("gen",1:length(genes_anotados_2)))))
colnames(tabla_de_abundancia) <- c("Abundancia","Gen")
tabla_de_abundancia
```

> 
  <p style="font-family: Comic Sans MS, serif; font-size:14pt">
  <span style="color:#004800">
  Por último (▼ abajo) podemos revisar si alguna lectura apareció en más de un gen. En este caso, no hubo ninguna secuencia que apareciera en más de 1, sin embargo hubo algunas lecturas que no aparecieron en ningún gen. Esto es debido a que estas lecturas probablemente corresponden a las lecturas finales. Si por ejemplo un gen va de la base 0 a 100, y otro 100 a 200, y una lectura va de las bases 94-108, entonces esa lectura no se registrará como parte de ningún gen. En datos reales de expresión de genes estas lecturas nunca se generarían (pues los RNA no se fusionan entre sí) pero es una artefacto pequeño.</p>
  
  ```{r import-data15, cache=TRUE, fig.width = 10, fig.asp = 1}
veces_apar <- NULL
for (i in 1:length(matriz_apariciones[,1])){
  veces_apar <- c(veces_apar,sum(matriz_apariciones[i,]))
}
any(veces_apar==2)
hist(veces_apar)

```
#####

###fastqc
###directorio root? hasta arriba/home/farell

###Kalisto
#conda activate
#kallisto comandos de classroom


#Sleuth



library("sleuth")

tx2gene <- function(){
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  return(t2g)
}

t2g <- tx2gene() #transcript to gene

base_dir <- "clase/26-apr/New folder/results" #DIrectorio donde vas a tener tus resultados, tomalo como si ./ fuera el r project

samples <- paste0("SRR4933", 66:71,"/kallisto") #el numbre del subdirectorio generado por kallisto con el tsv
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = rep(c("ctrl", "exp"), each=3), stringsAsFactors=FALSE) ###Generacion de grupos
so <- sleuth_prep(s2c, ~timepoint, target_mapping = t2g)
so <- sleuth_fit(so)
so <- sleuth_wt(so, which_beta="timepointt24h") 
sleuth_live(so)



#Binomial
dbinom(a,size=b,prob=c)
#a es el número de eventos "exitosos" (que salga cara)
#b es el número de eventos totales
#c es la probabilidad de que el evento a se suscite