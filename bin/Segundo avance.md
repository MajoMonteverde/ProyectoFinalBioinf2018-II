
# Streptanhus RAD-seq analyses

## by Majo Monteverde-Suárez

### Segundo avance de proyecto final del curso de Bioinfiormática 2018-2

Mis avances hasta el momento son:
1. Obtención del archivo de datos para analizar. Estos datos fueron generado por Nick Jensen (estudiante de alguna universidad en EUA que tiene contacto con mi tutora). Por el momento, sólo nos compartió los datos crudos y estoy esperando respuesta de correo para saber: ¿qué protocolo de RAD usó (RADoriginal, ddRAD, etc.)? ¿si son datos single o double paired RAD? y detalles sobre las enzimas utilizadas.
2. Estos datos fueron generados para 66 especies, que pertenecen a 26 géneros (aprox.) del complejo Thelypodieae (Brassicaceae)
3. Hasta el momento he empezado con el preprocesamiento de datos.
4. Demultiplexee y quité barcodes (falta verificar esto) las secuencias utilizando ipyrad.
5. Tengo los archivos de calidad generados por fastxtools (mean quality está entre 30 y 42)
6. He empezado a hacer el comando para correr las gráficas de calidad de fastxtools y verificar si realmente ya quité los barcodes y con base a esto, limpiar las secuencias.
7. Después de hacer las gráficas de calidad y probar al menos dos pipelines de filtrado de datos, pienso utilizar ipyrad para hacer algunos análisis filogenéticos. También pienso generar un archivo para poder hacer un análisis filogenético al menos en RAxML.


### Problemas a los que me he enfrentado
He tenido algunas complicaciones para poder usar las rutas en los forloops. Por ejemplo, en el paso 1.2.2 donde hago los archivos de calidad para las secuencias, en la parte del output no me dejaba guardarlos en la carpeta donde quería así que después los moví todos usando el comando mv.
Tengo ahorita algunas dudas sobre como correr el comando de las gráficas de calidad. Revisando el manual de FASTXtools el comando a usar es:
Usage: /usr/local/bin/fastq_quality_boxplot_graph.sh [-i INPUT.TXT] [-t TITLE] [-p] [-o OUTPUT]

así que tengo que ver bien sobre el script que según entiendo hay que utilizar y ver si me acepta el comando que puse como output para guardar los archivos o sucede algo similar ocmo lo que me pasó con los archivos de calidad.


### Otros avances:
He hecho una revisión profunda de métodos de RADseq, sus diferencias y visto artículos para comparar que han detallado sobre la limpieza de datos y los análisis que han hecho ¿pongo bibliografía de todo lo que he leído?
Tomé el curso de ipyrad impartido por Deren Eaton. El curso de bioinformática me ayudó un montón para poder entender bastante bien los comandos y explicación sobre ipyrad en el curso. Incluso creo que los que llevamos el curso sentimos mayor facilidad que las otras personas que no han tomado este curso.




### 1. Preprocessing

Take a look at the file in the .gz


```python
gunzip -c ./rad_data/1_preprocesamiento/UO_C601_1.fastq.gz | head -n 12
```


```python
@K00337:20:HF2H3BBXX:1:1101:1042:1191 1:N:0:1
NAACATGAAGTGCAGAGATTTGAGTGAAATTTGCTTTGGGTGAAAGTCCTCGTCCTTTACCTCTTGTGTTTTTACTGTTTGAGCGGTCTTCCCACCATTCTGGATATCCCACAACTTTGAAACATGAAGAGGCTTCATGTCCTGAGCGCCC
+
#AAFFJJJJFAFJJJJJFJJF<JFJJJJJJJJJJJAFJJJJFJFJJJJJFJJJJJJJJJJJJFJJJAJFJJJJJJJFJJJJJJJJF<FJFJJJJJJFJFJJJJJJF7JJJJFFFJFJAJAFJJJAAJFJJF<FFAAA<<JFJFJ7FFAJJJ
@K00337:20:HF2H3BBXX:1:1101:1062:1191 1:N:0:1
NTCCATAAGATGCAGACCACGATTAAGGCAGCGACGTCTCCGCGGGCGTATCAAAAGCCCGGGCTTAGGCCACCACCTTAATCCGCGTCGGTCCACGCCCCGAATCGATCGGCCGACCGCATCGCTCCGTTCCGCATCCGACCGGGACGCC
+
#--AAJFJAFAFJJJJJFJJFAFJFJJAAAJJAJJJJ<<A-7A7JJJJ<JJFJJFJJJFJJJJF-7A<FJJF-A777FJJF-FF-7AA7<JFJF<FF<7J<77<AF-A7-F7F-<<JJF))F)7<-7AF)<AAAF-AJJ)-<F<<))-)-)
@K00337:20:HF2H3BBXX:1:1101:1083:1191 1:N:0:1
NGTTGTCCGCTGCAGTGTACTCAGGTGCTTGTTAAGGGAATTGTCAAACTCGATAAGCTTGGTCGTGTAACCGGGTCTTTTGTAAAATTTCTTCCACTTTTTGATCTTCTGGTACTCGTCAACCATCTGCTCTCCCTTCTTCAGCAAATCC
+
#<<AFJJJJJAAJFJJJJJJFJJJAJJJJAFJJJJJF-<FJAF-FJJ<FJJFJJJJJJJFFJJJJJ<AFJJJAJA<A7F7FJJJJJFJF-FJJFJJJJJJ--77FFFJJJJJ<FF<J-FJ7<FFA-77<AAF<JJJJ-<FFJJFJJ7-A<F

```

### 1.1 Demultiplexing: remove barcodes

#### 1.1.1 Start ipyrad:
 -n create a new params file
 strep_demulti = name of the params file


```python
ipyrad -n strep_demulti

#  New file 'params-strep_demulti.txt' created in /home/majo/Escritorio/intento1_radstrep
```

#### 1.1.2 Visualize the params file


```python
cat params-strep_demulti.txt
```


```python
------- ipyrad params file (v.0.7.23)-------------------------------------------
strep_demulti                  ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
./                             ## [1] [project_dir]: Project dir (made in curdir if not present)
                               ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
                               ## [3] [barcodes_path]: Location of barcodes file
                               ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
denovo                         ## [5] [assembly_method]: Assembly method (denovo, reference, denovo+reference, denovo-reference)
                               ## [6] [reference_sequence]: Location of reference sequence file
rad                            ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
TGCAG,                         ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
33                             ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and very standard)
6                              ## [11] [mindepth_statistical]: Min depth for statistical base calling
6                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
10000                          ## [13] [maxdepth]: Max cluster depth within samples
0.85                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly
0                              ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
0                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
35                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
5, 5                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus (R1, R2)
8, 8                           ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus (R1, R2)
4                              ## [21] [min_samples_locus]: Min # samples per locus for output
20, 20                         ## [22] [max_SNPs_locus]: Max # SNPs per locus (R1, R2)
8, 8                           ## [23] [max_Indels_locus]: Max # of indels per locus (R1, R2)
0.5                            ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus (R1, R2)
0, 0, 0, 0                     ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)
0, 0, 0, 0                     ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)
p, s, v                        ## [27] [output_formats]: Output formats (see docs)
```

#### 1.1.3 Modify params file in step 2 and 3 which assign the path to the raw dara and to the pbarcodes file, respectively


```python
nano params-strep_demulti.txt
```


```python
 ------- ipyrad params file (v.0.7.23)-------------------------------------------
strep_demulti                  ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
./                             ## [1] [project_dir]: Project dir (made in curdir if not present)
./rad_data/1_preprocesamiento/UO_C601_1.fastq.gz   ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
./rada_data/barcodes/UO_C601_1_barcodes.txt        ## [3] [barcodes_path]: Location of barcodes file
```

### Step 1. Reads in barcodes file and the raw data. It separates each sample into a new independent file.


```python
## Run step 1, -p indicates the params file, -s indcates which step to run

ipyrad -p params-strep_demulti.txt -s 1
```


```python
 -------------------------------------------------------------
  ipyrad [v.0.7.23]
  Interactive assembly and analysis of RAD-seq data
 -------------------------------------------------------------
  New Assembly: strep_demulti
  establishing parallel connection:
  host compute node: [4 cores] on majo-Latitude-E7450
 
 Step 1: Demultiplexing fastq data to Samples
 [####################] 100%  chunking large files  | 1:03:00  
 [####################] 100%  sorting reads         | 1:16:50  
 [####################] 100%  writing/compressing   | 0:44:58
```


```python
## Informative results from currently executed steps

ipyrad -p params-strep_demulti.txt -r
```


```python
Summary stats of Assembly strep_demulti
------------------------------------------------
            state  reads_raw
DC_1            1      49168
DC_2            1     858987
DC_3            1     494565
DC_4            1     858275
FGXCONTROL      1     762703
ML_1            1    1828097
ML_2            1     327995
ML_3            1     190776
ML_4            1     131901
ML_5            1     145743
NJ_0003         1     979497
NJ_0005         1    1959192
NJ_0006         1    1014302
NJ_0007         1     966083
NJ_0008         1     984218
NJ_0009         1    2020934
NJ_0010         1    1955766
NJ_0011         1    1400084
NJ_0012         1     874498
NJ_0015         1     566369
NJ_0026         1    1070559
NJ_0029         1    2880290
NJ_0030         1    1619659
NJ_0031         1     257590
NJ_0032         1    1229165
NJ_0033         1    3024548
NJ_0034         1    3506655
NJ_0035         1     767502
NJ_0036         1    1354716
NJ_0037         1    1370969
NJ_0038         1    2941390
NJ_0039         1    2483993
NJ_0040         1    4031332
NJ_0041         1    5005922
NJ_0045         1    3326712
NJ_0047         1    7187677
NJ_0048         1     838577
NJ_0049         1    1820637
NJ_0051         1    1654916
NJ_0052         1    2807298
NJ_0053         1    4384643
NJ_0054         1    1844422
NJ_0055         1    2190132
NJ_0056         1    2799489
NJ_0057         1     745198
NJ_0058         1    3146868
NJ_0060         1     978737
NJ_0061         1    3770196
NJ_0062         1     159115
NJ_0063         1    3581430
NJ_0064         1    4964039
NJ_0065         1    1796649
NJ_0066         1     835816
NJ_0067         1     993597
NJ_0068         1    9412347
NJ_0069         1   22290334
NJ_0070         1   14961759
NJ_0071         1     889101
NJ_0072         1    2341715
NJ_0073         1    1433552
NJ_0074         1    7854099
NJ_0075         1    1009744
NJ_0076         1     882311
NJ_0077         1   14442857
NJ_0078         1    8557721
NJ_0079         1    1399754
NJ_0080         1     630119
NJ_0081         1    2318607
NJ_0082         1    8133790
NJ_0083         1    2131262
NJ_0084         1    2363286
NJ_0085         1    2475902
NJ_0086         1   11675566
NJ_0087         1    7009026
NJ_0088         1   12204090
NJ_0089         1    3957722
NM_3875         1    1113757
NM_759          1    2536858
NM_761          1     753500
NM_821          1    1289595
NM_839          1     632028
NM_863          1    1383160
NM_863D         1     976271
NM_864          1     998691
NM_864D         1     990931
NM_870          1    1131305
NM_871          1    1329276
NM_872          1    1041835
NM_873          1    1100105
NM_877          1     995419
NM_879          1     853565
NM_879D         1    2174081
NM_885          1    1147719
NM_SNBU         1     957172
NM_SNSA         1    2295401
NM_SNSA2        1     311295


Full stats files
------------------------------------------------
step 1: ./strep_demulti_fastqs/s1_demultiplex_stats.txt
step 2: None
step 3: None
step 4: None
step 5: None
step 6: None
step 7: None
```

### 1.2 Quality evaluation using Fastxtools

#### 1.2.1 Decompress all the files generated from Step 1.


```python
# Create a loop that decompress all the files in the strep_demulti_fastqs folder.

for i in strep_demulti_fastqs/*fastqc.gz; do
gunzip $i;
done
```

#### 1.2.2 Do FastX Statistics: fastx_quality_stats

Generate txt files with quality summary results


```python
for i in strep_demulti_fastqs/*.fastq; do
    fastx_quality_stats -i $i -o $i.txt;
    done
    
    
## Mean quality of samples is between 31 and 40
```

Create a directory to save output files


```python
mkdir strep_demulti_fastqs/1_fastx_quality_stats 
```

Move files of quality stats to the new folder


```python
mv strep_demulti_fastqs/*fastq.txt strep_demulti_fastqs/1_fastx_quality_stats
```

#### 1.2.3 Do fastq quality analysis


```python
for i in strep_demulti_fastqs/*.fastq; do
  fastqc $i;
    done
```

Create a directory to save the html files


```python
mkdir strep_demulti_fastqs/2_fastqc_qualitycharts
```

Move files of quality charts to the new folder


```python
mv strep_demulti_fastqs/*.html strep_demulti_fastqs/2_fastqc_qualitycharts
mv strep_demulti_fastqs/*.zip strep_demulti_fastqs/2_fastqc_qualitycharts
```

#### 1.2.4 Resume stats from fastqc quality charts (done in R version 3.4.4)


```python
############ Quality graph for Strptanthus sequences after demultiplexing ##############

## Load library
library(ggplot2)

## Read data
dat <- read.csv("../strep_demulti_fastqs/2_fastqc_qualitycharts/1_summary.csv")

## Plot
head (dat)
dat2<-table(dat$Wrong)
class(dat2)

## Remove row without data
dat2 <- as.data.frame(dat2)
dat2 <- dat2[-c(1),]

## Plot 
ggplot(data=dat2, aes(x=Var1, y= Freq)) + geom_bar(stat="identity")
myplot<-ggplot(data=dat2, aes(x=Var1, y= Freq))
myplot + geom_bar(stat="identity")

#Change colors
x <- ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) + geom_bar(stat="identity")
ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) +
  geom_bar(stat="identity", aes(color=Var1, fill=Var1))

## Remove grey background
y <- ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) +
  geom_bar(stat="identity", aes(color=Var1, fill=Var1)) + 
  theme_bw()

## Change axes names
z <- y + ylab("Frequency") + xlab("Problem")
z  

## Remove names for each bar
a <- z + theme(axis.text.x=element_blank())
a



sessionInfo()
R version 3.4.4 (2018-03-15)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.4 LTS

Matrix products: default
BLAS: /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=es_MX.UTF-8       LC_NUMERIC=C               LC_TIME=es_MX.UTF-8        LC_COLLATE=es_MX.UTF-8     LC_MONETARY=es_MX.UTF-8   
 [6] LC_MESSAGES=es_MX.UTF-8    LC_PAPER=es_MX.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=es_MX.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] rgeos_0.3-26             maptools_0.9-2           car_3.0-0                carData_3.0-1            plyr_1.8.4              
[6] ChemometricsWithR_0.1.11 ggplot2_2.2.1            sp_1.2-7                

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.16        magrittr_1.5        devtools_1.13.5     MASS_7.3-49         munsell_0.4.3       colorspace_1.3-2   
 [7] lattice_0.20-35     rlang_0.2.0         tools_3.4.4         grid_3.4.4          data.table_1.10.4-3 gtable_0.2.0       
[13] rio_0.5.10          withr_2.1.2         abind_1.4-5         readxl_1.1.0        yaml_2.1.18         lazyeval_0.2.1     
[19] digest_0.6.15       tibble_1.4.2        curl_3.2            memoise_1.1.0       haven_1.1.1         labeling_0.3       
[25] openxlsx_4.0.17     cellranger_1.1.0    compiler_3.4.4      pillar_1.2.2        forcats_0.3.0       pls_2.6-0          
[31] scales_0.5.0        foreign_0.8-69      kohonen_3.0.4      

```

# ![img](https://github.com/MajoMonteverde/ProyectoFinalBioinf2018-II/blob/master/figures/1_grafica_calidad.png)

Mayor problems found in the sequences after demultiplexing were:
1. Per tile sequence quality: reasons for seeing warnings or errors on this plot could be transient problems such as bubbles going through the flowcell, or they could be more permanent problems such as smudges on the flowcell or debris inside the flowcell lane.
2. Sequence duplication levels: high levels of duplication indicate some kind of enrichment bias (e.g. PCR over amplification)(FASTQC Manual)
3. Kmer content: any individually overrepresented sequences, even if not present at a high enough threshold to trigger the overrepresented sequences module will cause the Kmers from those sequences to be highly enriched in this module. 
4. Overrepresented sequences: this module will issue an error if any sequence is found to represent more than 1% of the total.

Based on this I decided to remove duplicated sequences for each individual.

### 1.3 Cleaning data

#### 1.3.1 To deal with overrepresented sequences, the fastq collapser option was used. This option collapses identical sequences in a FASTQ/A file into a single sequence (while maintaining reads counts). First, create a directory where collapsed sequences will be saved.


```python
mkdir strep_demulti_fastqs/3_collapsed
```


```python
## For one file
# fastx_collapser -v -i strep_demulti_fastqs/DC_1_R1_.fastq -o strep_demulti_fastqs/3_collapsed/collapsed_DC_1_R1 

for i in strep_demulti_fastqs/*.fastq; do
  fastx_collapser -v -i $i -o $i.txt;
done

```

Move collapsed files to 3_collapsed folder


```python
mv strep_demulti_fastqs/*.txt strep_demulti_fastqs/3_collapsed
```
