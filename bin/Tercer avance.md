
# Streptanhus RAD-seq analyses

## by Majo Monteverde-Suárez

### Tercer avance de proyecto final del curso de Bioinfiormática 2018-2

Mis avances hasta el momento son:
1. Terminé con la limpieza de datos: corrí fastqc.
2. He corrido hasta el paso 5 de 7 del tutorial de ipyrad utilizando mis datos genéticos.
3. He generado algunas gráficas en R para resumir algunos de los parámetros/resultados de la limpieza de datos
4. Empecé a organizar el repositorio. Hasta el momento tengo tres carpetas: bin donde he puesto los md de los avances y los actualicé en la lista de entrega. También dentro de esta carpeta hice una llamada scripts donde puse mis scripts de R. Las otras dos carpetas son data, aunque todavía no subo archivos ahí y figures donde puse las dos gráficas que he hecho hasta el momento.
5. [Liga al README](https://github.com/MajoMonteverde/ProyectoFinalBioinf2018-II/blob/master/README.md)


### Problemas a los que me he enfrentado
Los análisis llevan mucho tiempo en correr!!!!!


### Cosas que todavía planeo y/o me faltan por hacer
1. Quiero hacer algunas gráficas para resumir la info de los resultados de algunos pasos de ipyrad
2. Me falta hacer el corte de mi raw data a menos secuencias para poderlas subir a github. Aquí no estoy segura sobre que tanto de datos debo poner en el repositorio. Tengo muchos archivos muy pesados.
3. Generar al final un archivo que pueda correr en RAxML.
4. Para el análisis en RAxML, pienso correr un jmodeltest para decidir el modelo evolutivo y correr el análisis RAxML en Cypress gateway.

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
############ Quality graph for Streptanthus sequences after demultiplexing ##############

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


### 1.3 Cleaning data

#### 1.3.1 Remove low base sequence quality

Create a new folder to save trimmed sequences


```python
mkdir strep_demulti_fastqs/3_trimmed
```

For DC_1_R1, removed bases 84 to 136


```python
fastx_trimmer -f 1 -l 83 -i strep_demulti_fastqs/DC_1_R1_.fastq -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim
```

For all sequences remove 134 to 136


```python
for i in strep_demulti_fastqs/*.fastq; do
fastx_trimmer -f 1 -l 133 -i $i -o $i.txt;
done
```

Move trimmed files to 3_trimmed folder


```python
mv strep_demulti_fastqs/*q.txt strep_demulti_fastqs/3_trimmed/
```

For DC_3_R1, removed bases 120 to 133


```python
fastx_trimmer -f 1 -l 119 -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim.fastq -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim1.fastq
```

For ML_3_R1, removed bases 90 to 133


```python
fastx_trimmer -f 1 -l 89 -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim.fastq -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim1.fastq
```

For NM_839_R1, removed bases 120 to 136


```python
fastx_trimmer -f 1 -l 119 -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim.fastq -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim1.fastq
```

Verify fastq quality


```python
for i in strep_demulti_fastqs/3_trimmed/*.fastq; do
  fastqc $i;
    done
```

Create new folder to save fastq quality charts for trim files


```python
mkdir strep_demulti_fastqs/4_fastqc_trim_qualitycharts
```

Move files of quality charts to the new folder


```python
mv strep_demulti_fastqs/3_trimmed/*.html strep_demulti_fastqs/4_fastqc_trim_qualitycharts
mv strep_demulti_fastqs/3_trimmed/*.zip strep_demulti_fastqs/4_fastqc_trim_qualitycharts
```

#### 1.3.2 Resume stats from fastqc quality charts after trimming (done in R version 3.4.4)


```python
############ Quality graph for Streptanthus sequences after trimming ##############

## Load library
library(ggplot2)

## Read data
dat <- read.csv("../strep_demulti_fastqs/4_fastqc_trim_qualitycharts/1_summary_trim.csv")

## Plotting
head (dat)
dat2<-table(dat$Wrong)
class(dat2)

dat2 <- as.data.frame(dat2)
dat2 <- dat2[-c(1),]

## Plotting
ggplot(data=dat2, aes(x=Var1, y= Freq)) + geom_bar(stat="identity")
myplot<-ggplot(data=dat2, aes(x=Var1, y= Freq))
myplot + geom_bar(stat="identity")

#aes: change color
x <- ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) + geom_bar(stat="identity")
ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) +
  geom_bar(stat="identity", aes(color=Var1, fill=Var1))

##change grey background
y <- ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) +
  geom_bar(stat="identity", aes(color=Var1, fill=Var1)) + 
  theme_bw()

##change axes names
z <- y + ylab("Frequency") + xlab("Problem")
z  

##remove names from each bar
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
[1] ggplot2_2.2.1

loaded via a namespace (and not attached):
 [1] labeling_0.3     colorspace_1.3-2 scales_0.5.0     compiler_3.4.4   lazyeval_0.2.1   plyr_1.8.4       tools_3.4.4      pillar_1.2.2     gtable_0.2.0    
[10] tibble_1.4.2     yaml_2.1.19      Rcpp_0.12.16     grid_3.4.4       digest_0.6.15    rlang_0.2.0      munsell_0.4.3 
```

# ![img](https://github.com/MajoMonteverde/ProyectoFinalBioinf2018-II/blob/master/figures/1_grafica_calidad_trim.png)

Mayor problems found in the sequences after trimming were:
1. Per tile sequence quality
2. Sequence duplication levels
3. Kmer content
4. Overrepresented sequences

Also found some sequences with adapter content and with some low quality bases. Cleaned these data as follows:

#### 1.3.3 Remove adapter content (Fastx clipper) and clean low quality bases

For DC_1_R1, removed bases 84-133


```python
fastx_trimmer -f 1 -l 83 -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim.fastq -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_2
```

Remove adapter sequence (-a adapter sequence)


```python
Illumina Paired End PCR Primer 2

fastx_clipper -a CTGGTTCTGAGATAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_2 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_3 
fastx_clipper -a CTGGTTCTGAGATCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGAC -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_3 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_4 
fastx_clipper -a CTGGTTCTGAGATCGGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_4 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_5
fastx_clipper -a CTGGTTCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCG -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_5 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_6
fastx_clipper -a CTGGTTCTGAGAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGA -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_6 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_7
fastx_clipper -a CTGGTTCTGAGATCGGAAGAGAGATCGGAAGAGCGGTTCAGCAGGAATGC -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_7 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_8
fastx_clipper -a CTGGTTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGT -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_8 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_9
fastx_clipper -a CTGGTTCTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_9 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_10
fastx_clipper -a CTGGTTCTGAGATCGGAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGA -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_10 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_11
fastx_clipper -a CTGGTTCTGAGATCGGAAGAGATCGGAAGAGCGGTTCAGCAGGAATGCCG -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_11 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_12
fastx_clipper -a CTGGTTCTGAGATCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACC -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_12 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_13
fastx_clipper -a CTGGTTCTGAGATCGGAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_13 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_14
fastx_clipper -a CTGGTTCTGAGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGAT -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_14 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_15
fastx_clipper -a GAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAA -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_15 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_16
fastx_clipper -a CTGGTTCTGAGATCGGAAGAAGATCGGAAGAGCGGTTCAGCAGGAATGCC -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_16 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_17
fastx_clipper -a CTGGTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTA -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_17 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_18
fastx_clipper -a CTGGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTAT -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_18 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_19 
fastx_clipper -a CTGGTTCTGAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATC -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_19 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_20 
fastx_clipper -a CTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGC -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_20 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_21  
fastx_clipper -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCG -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_21 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_22
fastx_clipper -a CTGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATG -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_22 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_23
fastx_clipper -a CTGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_23 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_24
fastx_clipper -a CAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCC -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_24 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_25
fastx_clipper -a GAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAC -i strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_25 -o strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_26
```

For DC_3_R1, remove adapter sequence


```python
Illumina Paired End PCR Primer 2 
    
fastx_clipper -a ATTGGTCTTAGATCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGAC -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim1.fastq -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_2 
fastx_clipper -a ATTGGTCTTAGATAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_2 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_3 
fastx_clipper -a ATTGGTCTTAGATCGGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_3 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_4 
fastx_clipper -a ATTGGTCTTAGATCGGAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGA -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_4 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_5
fastx_clipper -a ATTGGTCTTAGATCGGAAGAGAGATCGGAAGAGCGGTTCAGCAGGAATGC -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_5 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_6
fastx_clipper -a ATTGGTCTTAGAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGA -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_6 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_7
fastx_clipper -a ATTGGTCTTAGATCGGAAGAGATCGGAAGAGCGGTTCAGCAGGAATGCCG -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_7 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_8
fastx_clipper -a ATTGGTCTTAGATCGGAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_8 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_9
fastx_clipper -a ATTGGTCTTAGATCGGAAGAAGATCGGAAGAGCGGTTCAGCAGGAATGCC -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_9 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_10
fastx_clipper -a ATTGGTCTTAGATCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACC -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_10 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_11 
fastx_clipper -a ATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTC -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_11 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_12 
fastx_clipper -a ATTGGTCTTAGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGAT -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_12 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_13
fastx_clipper -a ATTGGTCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCG -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_13 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_14
fastx_clipper -a ATTGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTC -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_14 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_15
fastx_clipper -a ATTGGTCTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_15 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_16
fastx_clipper -a ATTGGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTA -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_16 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_17
fastx_clipper -a ATTGGTCTTAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATC -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_17 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_18
fastx_clipper -a ATAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGC -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_18 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_19
fastx_clipper -a AAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCC -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_19 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_20
fastx_clipper -a ATTGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTAT -i strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_20 -o strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_21
```

For ML_3_R1, removed bases 86-89


```python
fastx_trimmer -f 1 -l 85 -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim1.fastq -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_2
```

Remove adapter sequence


```python
Illumina Paired End PCR Primer 2

fastx_clipper -a TAACGACCGAGATCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGAC -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_2 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_3
fastx_clipper -a TAACGACCGAGATAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_3 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_4
fastx_clipper -a TAACGACCGAGATCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACC -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_4 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_5
fastx_clipper -a TAACGACCGAGAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGA -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_5 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_6
fastx_clipper -a TAACGACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCG -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_6 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_7
fastx_clipper -a GAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAA -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_7 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_8
fastx_clipper -a TAACGACCGAGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGAT -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_8 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_9
fastx_clipper -a TAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGC -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_9 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_10
fastx_clipper -a TAACGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTA -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_10 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_11
fastx_clipper -a TAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATG -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_11 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_12
fastx_clipper -a TAACGACCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_12 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_13
fastx_clipper -a TAACGACCGAGATCGGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_13 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_14
fastx_clipper -a TAACGAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGT -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_14 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_15
fastx_clipper -a TAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCC -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_15 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_16
fastx_clipper -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCG -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_16 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_17
fastx_clipper -a TAACGACCGAGATCGGAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGA -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_17 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_18
fastx_clipper -a TAACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTAT -i strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_18 -o strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_19
```

For NJ_0075_R1, remove adapter sequence


```python
Illumina Paired End PCR Primer 2

fastx_clipper -a GTTAATGAGAGATCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGAC -i strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim.fastq -o strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_1
fastx_clipper -a GTTAATGAGAGATCGGAAGAGAGATCGGAAGAGCGGTTCAGCAGGAATGC -i strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_1 -o strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_2
fastx_clipper -a GTTAATGAGAGATAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG -i strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_2 -o strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_3
fastx_clipper -a GTTAATGAGAGAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGA -i strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_3 -o strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_4
fastx_clipper -a GTTAATGAGAGATCGGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA -i strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_4 -o strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_5
fastx_clipper -a GTTAATGAGAGATCGGAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGA -i strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_5 -o strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_6
fastx_clipper -a GTTAATGAGAGATCGGAAGAGATCGGAAGAGCGGTTCAGCAGGAATGCCG -i strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_6 -o strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_7
fastx_clipper -a GTTAATGAGAGATCGGAAGAAGATCGGAAGAGCGGTTCAGCAGGAATGCC -i strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_7 -o strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_8
fastx_clipper -a GTTAATGAGAGATCGGAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG -i strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_8 -o strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_9
fastx_clipper -a GTTAATGAGAGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGAT -i strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_9 -o strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_10
fastx_clipper -a GTTAATGAGAGATCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACC -i strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_10 -o strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_11
```

For NM_839_R1, removed bases 116-119


```python
fastx_trimmer -f 1 -l 85 -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim1.fastq -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_2
```

Remove adapter sequence


```python
Illumina Paired End PCR Primer 2

fastx_clipper -a GAGACGTGCAGATCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGAC -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_2 -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_3
fastx_clipper -a GAGACGTGCAGATCGGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_3 -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_4
fastx_clipper -a GAGACGTGCAGATCGGAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGA -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_4 -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_5
fastx_clipper -a GAGACGTGCAGATCGGAAGAGAGATCGGAAGAGCGGTTCAGCAGGAATGC -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_5 -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_6
fastx_clipper -a GAGACGTGCAGATAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_6 -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_7
fastx_clipper -a GAGACGTGCAGATCGGAAGAGATCGGAAGAGCGGTTCAGCAGGAATGCCG -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_7 -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_8
fastx_clipper -a GAGACGTGCAGATCGGAAGAAGATCGGAAGAGCGGTTCAGCAGGAATGCC -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_8 -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_9
fastx_clipper -a GAGACGTGCAGAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGA -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_9 -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_10
fastx_clipper -a GAGACGTGCAGATCGGAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_10 -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_11
fastx_clipper -a GAGACGTGCAGATCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACC -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_11 -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_12
fastx_clipper -a GAGACGTGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_12 -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_13
fastx_clipper -a GAGACGTGCAGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGAT -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_13 -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_14
fastx_clipper -a GAGACGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGT -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_14 -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_15
fastx_clipper -a GAGACGTGCAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATC -i strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_15 -o strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_16
```

For NM_864_R1, remove adapter sequence


```python
Illumina Paired End PCR Primer 2

fastx_clipper -a CGGACAACAAGATCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGAC -i strep_demulti_fastqs/3_trimmed/NM_864_R1_trim.fastq -o strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_1
fastx_clipper -a CGGACAACAAGATCGGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA -i strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_1 -o strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_2
fastx_clipper -a CGGACAACAAGATCGGAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGA -i strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_2 -o strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_3
fastx_clipper -a CGGACAACAAGATCGGAAGAGAGATCGGAAGAGCGGTTCAGCAGGAATGC -i strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_3 -o strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_4
fastx_clipper -a CGGACAACAAGATAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG -i strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_4 -o strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_5
fastx_clipper -a CGGACAACAAGAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGA -i strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_5 -o strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_6
fastx_clipper -a CGGACAACAAGATCGGAAGAAGATCGGAAGAGCGGTTCAGCAGGAATGCC -i strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_6 -o strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_7
fastx_clipper -a CGGACAACAAGATCGGAAGAGATCGGAAGAGCGGTTCAGCAGGAATGCCG -i strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_7 -o strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_8
fastx_clipper -a CGGACAACAAGATCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACC -i strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_8 -o strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_9
fastx_clipper -a CGGACAACAAGATCGGAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG -i strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_9 -o strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_10
fastx_clipper -a CGGACAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCG -i strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_10 -o strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_11
fastx_clipper -a CGGACAACAAGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGAT -i strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_11 -o strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_12
fastx_clipper -a CGGACAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGT -i strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_12 -o strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_13
fastx_clipper -a CGGACAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTA -i strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_13 -o strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_14
```

Do fastqc analysis for DC_1_R1_trim_26, DC_3_R1_trim_21, ML_3_R1_trim_19, NJ_0075_R1_trim_11, NM_839_R1_trim_16, NM_864_R1_trim_14


```python
fastqc strep_demulti_fastqs/3_trimmed/DC_1_R1_trim_26
fastqc strep_demulti_fastqs/3_trimmed/DC_3_R1_trim_21
fastqc strep_demulti_fastqs/3_trimmed/ML_3_R1_trim_19
fastqc strep_demulti_fastqs/3_trimmed/NJ_0075_R1_trim_11
fastqc strep_demulti_fastqs/3_trimmed/NM_839_R1_trim_16
fastqc strep_demulti_fastqs/3_trimmed/NM_864_R1_trim_14
```

Create new directory for new fastqc results


```python
mkdir strep_demulti_fastqs/5_fastqc_clipper
```

Move fastqc result files to strep_demulti_fastqs/4_fastqc_clipper


```python
mv strep_demulti_fastqs/3_trimmed/*.html strep_demulti_fastqs/5_fastqc_clipper
mv strep_demulti_fastqs/3_trimmed/*.zip strep_demulti_fastqs/5_fastqc_clipper
```

Adapter content is good, length sequence with warning.

### 2. Filtering data in ipyrad based on sequence quality

#### 2.1 Filter reads using ipyrad

Change params files


```python
nano params-strep_demulti.txt
```


```python
strep_demulti                  ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
./                             ## [1] [project_dir]: Project dir (made in curdir if not present)
./rad_data/1_preprocesamiento/UO_C601_1.fastq.gz   ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
./rad_data/Barcodes/UO_C601_1_barcodes.txt         ## [3] [barcodes_path]: Location of barcodes file
../intento1_radstrep/strep_demulti_fastqs/*.gz     ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
denovo                         ## [5] [assembly_method]: Assembly method (denovo, reference, denovo+reference, denovo-reference)

```

Run step 2 of ipyrad


```python
ipyrad -p params-strep_demulti.txt -s 2
```

Current stats after filtering


```python
ipyrad -p params-strep_demulti.txt -r
```


```python
Summary stats of Assembly strep_demulti
------------------------------------------------
            state  reads_raw  reads_passed_filter
DC_1            2      49168                 4248
DC_2            2     858987               856604
DC_3            2     494565               164572
DC_4            2     858275               855749
FGXCONTROL      2     762703               761073
ML_1            2    1828097              1820777
ML_2            2     327995               326565
ML_3            2     190776                48843
ML_4            2     131901               131054
ML_5            2     145743               144668
NJ_0003         2     979497               957742
NJ_0005         2    1959192              1938568
NJ_0006         2    1014302              1011300
NJ_0007         2     966083               963291
NJ_0008         2     984218               981660
NJ_0009         2    2020934              2014858
NJ_0010         2    1955766              1937965
NJ_0011         2    1400084              1396528
NJ_0012         2     874498               871261
NJ_0015         2     566369               565033
NJ_0026         2    1070559              1067814
NJ_0029         2    2880290              2870481
NJ_0030         2    1619659              1602807
NJ_0031         2     257590               256514
NJ_0032         2    1229165              1212594
NJ_0033         2    3024548              3016220
NJ_0034         2    3506655              3496989
NJ_0035         2     767502               765200
NJ_0036         2    1354716              1350229
NJ_0037         2    1370969              1366922
NJ_0038         2    2941390              2930885
NJ_0039         2    2483993              2476015
NJ_0040         2    4031332              4000726
NJ_0041         2    5005922              4993123
NJ_0045         2    3326712              3317427
NJ_0047         2    7187677              7163615
NJ_0048         2     838577               834492
NJ_0049         2    1820637              1783485
NJ_0051         2    1654916              1648948
NJ_0052         2    2807298              2797181
NJ_0053         2    4384643              4373856
NJ_0054         2    1844422              1839502
NJ_0055         2    2190132              2183890
NJ_0056         2    2799489              2791712
NJ_0057         2     745198               727419
NJ_0058         2    3146868              3125059
NJ_0060         2     978737               976257
NJ_0061         2    3770196              3760323
NJ_0062         2     159115               158636
NJ_0063         2    3581430              3572356
NJ_0064         2    4964039              4951402
NJ_0065         2    1796649              1773103
NJ_0066         2     835816               833685
NJ_0067         2     993597               974608
NJ_0068         2    9412347              9388139
NJ_0069         2   22290334             22236146
NJ_0070         2   14961759             14926921
NJ_0071         2     889101               887109
NJ_0072         2    2341715              2335426
NJ_0073         2    1433552              1415310
NJ_0074         2    7854099              7833293
NJ_0075         2    1009744               761190
NJ_0076         2     882311               846632
NJ_0077         2   14442857             14410267
NJ_0078         2    8557721              8535933
NJ_0079         2    1399754              1396580
NJ_0080         2     630119               628747
NJ_0081         2    2318607              2294005
NJ_0082         2    8133790              8097476
NJ_0083         2    2131262              2126143
NJ_0084         2    2363286              2356677
NJ_0085         2    2475902              2470464
NJ_0086         2   11675566             11651298
NJ_0087         2    7009026              6993284
NJ_0088         2   12204090             12165858
NJ_0089         2    3957722              3950461
NM_3875         2    1113757              1110586
NM_759          2    2536858              2501458
NM_761          2     753500               753323
NM_821          2    1289595              1289300
NM_839          2     632028               398312
NM_863          2    1383160              1378601
NM_863D         2     976271               976084
NM_864          2     998691               672885
NM_864D         2     990931               990721
NM_870          2    1131305              1131074
NM_871          2    1329276              1325158
NM_872          2    1041835              1041619
NM_873          2    1100105              1099872
NM_877          2     995419               992256
NM_879          2     853565               853389
NM_879D         2    2174081              2173587
NM_885          2    1147719              1144390
NM_SNBU         2     957172               956963
NM_SNSA         2    2295401              2294944
NM_SNSA2        2     311295               311229


Full stats files
------------------------------------------------
step 1: ./strep_demulti_fastqs/s1_demultiplex_stats.txt
step 2: ./strep_demulti_edits/s2_rawedit_stats.txt
step 3: None
step 4: None
step 5: None
step 6: None
step 7: None
```

# Aquí pienso hacer una gráfica-resumen como de porcentajes/promedio de cuántas lecturas pasaron por el filtro

### 3. Clustering data

This step de-replicates and clusters reads within each sample with 0.85 clustering thershold 


```python
ipyrad -p params-strep_demulti.txt -s 3
```


```python
/home/majo/miniconda2/lib/python2.7/site-packages/h5py/__init__.py:36: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.
  from ._conv import register_converters as _register_converters

 -------------------------------------------------------------
  ipyrad [v.0.7.23]
  Interactive assembly and analysis of RAD-seq data
 -------------------------------------------------------------
  loading Assembly: strep_demulti
  from saved path: ~/Escritorio/intento1_radstrep/strep_demulti.json
  establishing parallel connection:
  host compute node: [4 cores] on majo-Latitude-E7450

  Step 3: Clustering/Mapping reads
  [####################] 100%  dereplicating         | 4:43:22  
  [####################] 100%  clustering            | 4:17:02  
  [####################] 100%  building clusters     | 0:00:21  
  [####################] 100%  chunking              | 0:00:02  
  [####################] 100%  aligning              | 9:50:12  
  [####################] 100%  concatenating         | 0:00:34 
```


```python
ipyrad -p params-strep_demulti.txt -r
```


```python
/home/majo/miniconda2/lib/python2.7/site-packages/h5py/__init__.py:36: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.
  from ._conv import register_converters as _register_converters

Summary stats of Assembly strep_demulti
------------------------------------------------
            state  reads_raw  reads_passed_filter  clusters_total  \
DC_1            3      49168                 4248             983   
DC_2            3     858987               856604           81521   
DC_3            3     494565               164572           21655   
DC_4            3     858275               855749           67836   
FGXCONTROL      3     762703               761073            7306   
ML_1            3    1828097              1820777          147738   
ML_2            3     327995               326565           79782   
ML_3            3     190776                48843           13427   
ML_4            3     131901               131054           33471   
ML_5            3     145743               144668           35598   
NJ_0003         3     979497               957742           72071   
NJ_0005         3    1959192              1938568           71159   
NJ_0006         3    1014302              1011300           67406   
NJ_0007         3     966083               963291           75635   
NJ_0008         3     984218               981660           69616   
NJ_0009         3    2020934              2014858          106717   
NJ_0010         3    1955766              1937965          101886   
NJ_0011         3    1400084              1396528           65129   
NJ_0012         3     874498               871261           71150   
NJ_0015         3     566369               565033           59009   
NJ_0026         3    1070559              1067814           80653   
NJ_0029         3    2880290              2870481          154090   
NJ_0030         3    1619659              1602807           78155   
NJ_0031         3     257590               256514           37702   
NJ_0032         3    1229165              1212594           66432   
NJ_0033         3    3024548              3016220           79836   
NJ_0034         3    3506655              3496989          101069   
NJ_0035         3     767502               765200           63764   
NJ_0036         3    1354716              1350229           78747   
NJ_0037         3    1370969              1366922           84087   
NJ_0038         3    2941390              2930885           84697   
NJ_0039         3    2483993              2476015           86744   
NJ_0040         3    4031332              4000726           94680   
NJ_0041         3    5005922              4993123           93988   
NJ_0045         3    3326712              3317427           99826   
NJ_0047         3    7187677              7163615           98509   
NJ_0048         3     838577               834492           61587   
NJ_0049         3    1820637              1783485           75416   
NJ_0051         3    1654916              1648948           78887   
NJ_0052         3    2807298              2797181           78465   
NJ_0053         3    4384643              4373856           93772   
NJ_0054         3    1844422              1839502           69405   
NJ_0055         3    2190132              2183890          106435   
NJ_0056         3    2799489              2791712          104933   
NJ_0057         3     745198               727419           68637   
NJ_0058         3    3146868              3125059          103023   
NJ_0060         3     978737               976257           56523   
NJ_0061         3    3770196              3760323           75802   
NJ_0062         3     159115               158636           28172   
NJ_0063         3    3581430              3572356          178327   
NJ_0064         3    4964039              4951402           91954   
NJ_0065         3    1796649              1773103           50831   
NJ_0066         3     835816               833685           61751   
NJ_0067         3     993597               974608           77524   
NJ_0068         3    9412347              9388139          187474   
NJ_0069         3   22290334             22236146          248468   
NJ_0070         3   14961759             14926921          301164   
NJ_0071         3     889101               887109           68052   
NJ_0072         3    2341715              2335426           73032   
NJ_0073         3    1433552              1415310           84120   
NJ_0074         3    7854099              7833293          128223   
NJ_0075         3    1009744               761190           56479   
NJ_0076         3     882311               846632           76321   
NJ_0077         3   14442857             14410267          116013   
NJ_0078         3    8557721              8535933           95483   
NJ_0079         3    1399754              1396580          138284   
NJ_0080         3     630119               628747           69109   
NJ_0081         3    2318607              2294005          128507   
NJ_0082         3    8133790              8097476          173446   
NJ_0083         3    2131262              2126143           84477   
NJ_0084         3    2363286              2356677           86809   
NJ_0085         3    2475902              2470464          112703   
NJ_0086         3   11675566             11651298          139996   
NJ_0087         3    7009026              6993284          100563   
NJ_0088         3   12204090             12165858          154044   
NJ_0089         3    3957722              3950461          112799   
NM_3875         3    1113757              1110586           56976   
NM_759          3    2536858              2501458          115954   
NM_761          3     753500               753323           89746   
NM_821          3    1289595              1289300           99947   
NM_839          3     632028               398312           41966   
NM_863          3    1383160              1378601           83393   
NM_863D         3     976271               976084           89355   
NM_864          3     998691               672885           97158   
NM_864D         3     990931               990721          139863   
NM_870          3    1131305              1131074           62898   
NM_871          3    1329276              1325158           68275   
NM_872          3    1041835              1041619           89021   
NM_873          3    1100105              1099872           64867   
NM_877          3     995419               992256           57088   
NM_879          3     853565               853389           56418   
NM_879D         3    2174081              2173587           68398   
NM_885          3    1147719              1144390           57111   
NM_SNBU         3     957172               956963           61917   
NM_SNSA         3    2295401              2294944          234173   
NM_SNSA2        3     311295               311229           63965   

            clusters_hidepth  
DC_1                     134  
DC_2                   40564  
DC_3                    8940  
DC_4                   36749  
FGXCONTROL              4728  
ML_1                   80117  
ML_2                   18778  
ML_3                    2348  
ML_4                    6493  
ML_5                    7615  
NJ_0003                42230  
NJ_0005                54567  
NJ_0006                41121  
NJ_0007                42620  
NJ_0008                41098  
NJ_0009                66563  
NJ_0010                63185  
NJ_0011                47935  
NJ_0012                40850  
NJ_0015                28088  
NJ_0026                45043  
NJ_0029                75645  
NJ_0030                50793  
NJ_0031                13496  
NJ_0032                45314  
NJ_0033                59051  
NJ_0034                76757  
NJ_0035                34729  
NJ_0036                51253  
NJ_0037                49383  
NJ_0038                60158  
NJ_0039                62243  
NJ_0040                70184  
NJ_0041                73180  
NJ_0045                70981  
NJ_0047                70956  
NJ_0048                39186  
NJ_0049                54130  
NJ_0051                52604  
NJ_0052                58132  
NJ_0053                72883  
NJ_0054                53049  
NJ_0055                67244  
NJ_0056                72063  
NJ_0057                32991  
NJ_0058                69527  
NJ_0060                38635  
NJ_0061                59500  
NJ_0062                 5953  
NJ_0063                86210  
NJ_0064                67246  
NJ_0065                18537  
NJ_0066                32084  
NJ_0067                30886  
NJ_0068               104428  
NJ_0069               108354  
NJ_0070               136434  
NJ_0071                37385  
NJ_0072                55911  
NJ_0073                50145  
NJ_0074                88029  
NJ_0075                29780  
NJ_0076                38313  
NJ_0077                77325  
NJ_0078                71618  
NJ_0079                60260  
NJ_0080                24387  
NJ_0081                72663  
NJ_0082               118211  
NJ_0083                60958  
NJ_0084                64307  
NJ_0085                62525  
NJ_0086                87387  
NJ_0087                72859  
NJ_0088                93824  
NJ_0089                45170  
NM_3875                37545  
NM_759                 57548  
NM_761                 37086  
NM_821                 48760  
NM_839                 21401  
NM_863                 51115  
NM_863D                42406  
NM_864                 35553  
NM_864D                49943  
NM_870                 40324  
NM_871                 43697  
NM_872                 46563  
NM_873                 41092  
NM_877                 37448  
NM_879                 34371  
NM_879D                46517  
NM_885                 37250  
NM_SNBU                35248  
NM_SNSA               117190  
NM_SNSA2               17634  


Full stats files
------------------------------------------------
step 1: ./strep_demulti_fastqs/s1_demultiplex_stats.txt
step 2: ./strep_demulti_edits/s2_rawedit_stats.txt
step 3: ./strep_demulti_clust_0.85/s3_cluster_stats.txt
step 4: None
step 5: None
step 6: None
step 7: None

```

### 4. Joint estimation of heterozygosity and error rate

Jointly estimates sequencing error and heterozygosity to disentangle which reads are "real" and which are sequencing error. 


```python
ipyrad -p params-strep_demulti.txt -s 4
```


```python
/home/majo/miniconda2/lib/python2.7/site-packages/h5py/__init__.py:36: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.
  from ._conv import register_converters as _register_converters

 -------------------------------------------------------------
  ipyrad [v.0.7.23]
  Interactive assembly and analysis of RAD-seq data
 -------------------------------------------------------------
  loading Assembly: strep_demulti
  from saved path: ~/Escritorio/intento1_radstrep/strep_demulti.json
  establishing parallel connection:
  host compute node: [4 cores] on majo-Latitude-E7450

  Step 4: Joint estimation of error rate and heterozygosity
  [####################] 100%  inferring [H, E]      | 1:02:34  
```


```python
ipyrad -p params-strep_demulti.txt -r
```


```python
/home/majo/miniconda2/lib/python2.7/site-packages/h5py/__init__.py:36: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.
  from ._conv import register_converters as _register_converters

Summary stats of Assembly strep_demulti
------------------------------------------------
            state  reads_raw  reads_passed_filter  clusters_total  \
DC_1            4      49168                 4248             983   
DC_2            4     858987               856604           81521   
DC_3            4     494565               164572           21655   
DC_4            4     858275               855749           67836   
FGXCONTROL      4     762703               761073            7306   
ML_1            4    1828097              1820777          147738   
ML_2            4     327995               326565           79782   
ML_3            4     190776                48843           13427   
ML_4            4     131901               131054           33471   
ML_5            4     145743               144668           35598   
NJ_0003         4     979497               957742           72071   
NJ_0005         4    1959192              1938568           71159   
NJ_0006         4    1014302              1011300           67406   
NJ_0007         4     966083               963291           75635   
NJ_0008         4     984218               981660           69616   
NJ_0009         4    2020934              2014858          106717   
NJ_0010         4    1955766              1937965          101886   
NJ_0011         4    1400084              1396528           65129   
NJ_0012         4     874498               871261           71150   
NJ_0015         4     566369               565033           59009   
NJ_0026         4    1070559              1067814           80653   
NJ_0029         4    2880290              2870481          154090   
NJ_0030         4    1619659              1602807           78155   
NJ_0031         4     257590               256514           37702   
NJ_0032         4    1229165              1212594           66432   
NJ_0033         4    3024548              3016220           79836   
NJ_0034         4    3506655              3496989          101069   
NJ_0035         4     767502               765200           63764   
NJ_0036         4    1354716              1350229           78747   
NJ_0037         4    1370969              1366922           84087   
NJ_0038         4    2941390              2930885           84697   
NJ_0039         4    2483993              2476015           86744   
NJ_0040         4    4031332              4000726           94680   
NJ_0041         4    5005922              4993123           93988   
NJ_0045         4    3326712              3317427           99826   
NJ_0047         4    7187677              7163615           98509   
NJ_0048         4     838577               834492           61587   
NJ_0049         4    1820637              1783485           75416   
NJ_0051         4    1654916              1648948           78887   
NJ_0052         4    2807298              2797181           78465   
NJ_0053         4    4384643              4373856           93772   
NJ_0054         4    1844422              1839502           69405   
NJ_0055         4    2190132              2183890          106435   
NJ_0056         4    2799489              2791712          104933   
NJ_0057         4     745198               727419           68637   
NJ_0058         4    3146868              3125059          103023   
NJ_0060         4     978737               976257           56523   
NJ_0061         4    3770196              3760323           75802   
NJ_0062         4     159115               158636           28172   
NJ_0063         4    3581430              3572356          178327   
NJ_0064         4    4964039              4951402           91954   
NJ_0065         4    1796649              1773103           50831   
NJ_0066         4     835816               833685           61751   
NJ_0067         4     993597               974608           77524   
NJ_0068         4    9412347              9388139          187474   
NJ_0069         4   22290334             22236146          248468   
NJ_0070         4   14961759             14926921          301164   
NJ_0071         4     889101               887109           68052   
NJ_0072         4    2341715              2335426           73032   
NJ_0073         4    1433552              1415310           84120   
NJ_0074         4    7854099              7833293          128223   
NJ_0075         4    1009744               761190           56479   
NJ_0076         4     882311               846632           76321   
NJ_0077         4   14442857             14410267          116013   
NJ_0078         4    8557721              8535933           95483   
NJ_0079         4    1399754              1396580          138284   
NJ_0080         4     630119               628747           69109   
NJ_0081         4    2318607              2294005          128507   
NJ_0082         4    8133790              8097476          173446   
NJ_0083         4    2131262              2126143           84477   
NJ_0084         4    2363286              2356677           86809   
NJ_0085         4    2475902              2470464          112703   
NJ_0086         4   11675566             11651298          139996   
NJ_0087         4    7009026              6993284          100563   
NJ_0088         4   12204090             12165858          154044   
NJ_0089         4    3957722              3950461          112799   
NM_3875         4    1113757              1110586           56976   
NM_759          4    2536858              2501458          115954   
NM_761          4     753500               753323           89746   
NM_821          4    1289595              1289300           99947   
NM_839          4     632028               398312           41966   
NM_863          4    1383160              1378601           83393   
NM_863D         4     976271               976084           89355   
NM_864          4     998691               672885           97158   
NM_864D         4     990931               990721          139863   
NM_870          4    1131305              1131074           62898   
NM_871          4    1329276              1325158           68275   
NM_872          4    1041835              1041619           89021   
NM_873          4    1100105              1099872           64867   
NM_877          4     995419               992256           57088   
NM_879          4     853565               853389           56418   
NM_879D         4    2174081              2173587           68398   
NM_885          4    1147719              1144390           57111   
NM_SNBU         4     957172               956963           61917   
NM_SNSA         4    2295401              2294944          234173   
NM_SNSA2        4     311295               311229           63965   

            clusters_hidepth  hetero_est  error_est  
DC_1                     134    0.005308   0.001329  
DC_2                   40564    0.016693   0.002049  
DC_3                    8940    0.005010   0.000918  
DC_4                   36749    0.007573   0.001608  
FGXCONTROL              4728    0.001282   0.001397  
ML_1                   80117    0.013519   0.002079  
ML_2                   18778    0.009793   0.001957  
ML_3                    2348    0.009088   0.001667  
ML_4                    6493    0.012995   0.002980  
ML_5                    7615    0.011525   0.002637  
NJ_0003                42230    0.021039   0.003053  
NJ_0005                54567    0.016550   0.001732  
NJ_0006                41121    0.017979   0.001936  
NJ_0007                42620    0.025856   0.002690  
NJ_0008                41098    0.019840   0.002787  
NJ_0009                66563    0.030777   0.003647  
NJ_0010                63185    0.029828   0.003458  
NJ_0011                47935    0.017666   0.001842  
NJ_0012                40850    0.020030   0.002992  
NJ_0015                28088    0.020467   0.002758  
NJ_0026                45043    0.022008   0.002658  
NJ_0029                75645    0.018572   0.001948  
NJ_0030                50793    0.021185   0.002645  
NJ_0031                13496    0.018112   0.002326  
NJ_0032                45314    0.019884   0.001974  
NJ_0033                59051    0.019101   0.001759  
NJ_0034                76757    0.019749   0.002993  
NJ_0035                34729    0.019773   0.002375  
NJ_0036                51253    0.022772   0.002636  
NJ_0037                49383    0.021730   0.002497  
NJ_0038                60158    0.021085   0.001795  
NJ_0039                62243    0.020848   0.002412  
NJ_0040                70184    0.017455   0.001689  
NJ_0041                73180    0.023126   0.002051  
NJ_0045                70981    0.024539   0.002568  
NJ_0047                70956    0.021202   0.002121  
NJ_0048                39186    0.019895   0.001952  
NJ_0049                54130    0.019663   0.002080  
NJ_0051                52604    0.022210   0.002460  
NJ_0052                58132    0.017298   0.001821  
NJ_0053                72883    0.020625   0.002306  
NJ_0054                53049    0.018328   0.001846  
NJ_0055                67244    0.030217   0.003633  
NJ_0056                72063    0.027690   0.002763  
NJ_0057                32991    0.020847   0.002424  
NJ_0058                69527    0.022007   0.002589  
NJ_0060                38635    0.008429   0.002897  
NJ_0061                59500    0.018336   0.001663  
NJ_0062                 5953    0.015460   0.002118  
NJ_0063                86210    0.014303   0.001939  
NJ_0064                67246    0.022482   0.001792  
NJ_0065                18537    0.021576   0.005675  
NJ_0066                32084    0.020127   0.003021  
NJ_0067                30886    0.027449   0.007137  
NJ_0068               104428    0.019692   0.002504  
NJ_0069               108354    0.018701   0.001270  
NJ_0070               136434    0.017976   0.001526  
NJ_0071                37385    0.022200   0.002734  
NJ_0072                55911    0.020707   0.001863  
NJ_0073                50145    0.026550   0.002836  
NJ_0074                88029    0.027581   0.002476  
NJ_0075                29780    0.018548   0.002437  
NJ_0076                38313    0.022905   0.003287  
NJ_0077                77325    0.018676   0.001548  
NJ_0078                71618    0.017639   0.001730  
NJ_0079                60260    0.014759   0.002580  
NJ_0080                24387    0.027919   0.008476  
NJ_0081                72663    0.028462   0.005226  
NJ_0082               118211    0.032335   0.005040  
NJ_0083                60958    0.019438   0.003130  
NJ_0084                64307    0.021354   0.002587  
NJ_0085                62525    0.019904   0.002479  
NJ_0086                87387    0.025240   0.002026  
NJ_0087                72859    0.021598   0.002043  
NJ_0088                93824    0.021661   0.002133  
NJ_0089                45170    0.018329   0.003797  
NM_3875                37545    0.008567   0.001674  
NM_759                 57548    0.012299   0.002059  
NM_761                 37086    0.009551   0.002154  
NM_821                 48760    0.011140   0.002314  
NM_839                 21401    0.009090   0.001479  
NM_863                 51115    0.009782   0.002075  
NM_863D                42406    0.011604   0.002247  
NM_864                 35553    0.008550   0.001873  
NM_864D                49943    0.011054   0.002483  
NM_870                 40324    0.009842   0.002176  
NM_871                 43697    0.009610   0.001973  
NM_872                 46563    0.008602   0.002237  
NM_873                 41092    0.011113   0.001980  
NM_877                 37448    0.009946   0.001902  
NM_879                 34371    0.009148   0.002134  
NM_879D                46517    0.010915   0.002089  
NM_885                 37250    0.008439   0.001619  
NM_SNBU                35248    0.009262   0.002337  
NM_SNSA               117190    0.018197   0.002748  
NM_SNSA2               17634    0.007853   0.002339  


Full stats files
------------------------------------------------
step 1: ./strep_demulti_fastqs/s1_demultiplex_stats.txt
step 2: ./strep_demulti_edits/s2_rawedit_stats.txt
step 3: ./strep_demulti_clust_0.85/s3_cluster_stats.txt
step 4: ./strep_demulti_clust_0.85/s4_joint_estimate.txt
step 5: None
step 6: None
step 7: None

```

### 5. Consensus base calls

Uses de inferred error rate and heterozygosity to call the consesnsus of sequences within each cluster: identify real haplotypes at each locus within each sample


```python
ipyrad -p params-strep_demulti.txt -s 5
```


```python
ipyrad -p params-strep_demulti.txt -s 5
/home/majo/miniconda2/lib/python2.7/site-packages/h5py/__init__.py:36: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.
  from ._conv import register_converters as _register_converters

 -------------------------------------------------------------
  ipyrad [v.0.7.23]
  Interactive assembly and analysis of RAD-seq data
 -------------------------------------------------------------
  loading Assembly: strep_demulti
  from saved path: ~/Escritorio/intento1_radstrep/strep_demulti.json
  establishing parallel connection:
  host compute node: [4 cores] on majo-Latitude-E7450

  Step 5: Consensus base calling 
  Mean error  [0.00248 sd=0.00109]
  Mean hetero [0.01768 sd=0.00656]
  [####################] 100%  calculating depths    | 0:02:14  
  [####################] 100%  chunking clusters     | 0:03:11  
  [####################] 100%  consens calling       | 3:51:37
```


```python
ipyrad -p params-strep_demulti.txt -r
```


```python
ipyrad -p params-strep_demulti.txt -r
/home/majo/miniconda2/lib/python2.7/site-packages/h5py/__init__.py:36: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.
  from ._conv import register_converters as _register_converters

Summary stats of Assembly strep_demulti
------------------------------------------------
            state  reads_raw  reads_passed_filter  clusters_total  \
DC_1            5      49168                 4248             983   
DC_2            5     858987               856604           81521   
DC_3            5     494565               164572           21655   
DC_4            5     858275               855749           67836   
FGXCONTROL      5     762703               761073            7306   
ML_1            5    1828097              1820777          147738   
ML_2            5     327995               326565           79782   
ML_3            5     190776                48843           13427   
ML_4            5     131901               131054           33471   
ML_5            5     145743               144668           35598   
NJ_0003         5     979497               957742           72071   
NJ_0005         5    1959192              1938568           71159   
NJ_0006         5    1014302              1011300           67406   
NJ_0007         5     966083               963291           75635   
NJ_0008         5     984218               981660           69616   
NJ_0009         5    2020934              2014858          106717   
NJ_0010         5    1955766              1937965          101886   
NJ_0011         5    1400084              1396528           65129   
NJ_0012         5     874498               871261           71150   
NJ_0015         5     566369               565033           59009   
NJ_0026         5    1070559              1067814           80653   
NJ_0029         5    2880290              2870481          154090   
NJ_0030         5    1619659              1602807           78155   
NJ_0031         5     257590               256514           37702   
NJ_0032         5    1229165              1212594           66432   
NJ_0033         5    3024548              3016220           79836   
NJ_0034         5    3506655              3496989          101069   
NJ_0035         5     767502               765200           63764   
NJ_0036         5    1354716              1350229           78747   
NJ_0037         5    1370969              1366922           84087   
NJ_0038         5    2941390              2930885           84697   
NJ_0039         5    2483993              2476015           86744   
NJ_0040         5    4031332              4000726           94680   
NJ_0041         5    5005922              4993123           93988   
NJ_0045         5    3326712              3317427           99826   
NJ_0047         5    7187677              7163615           98509   
NJ_0048         5     838577               834492           61587   
NJ_0049         5    1820637              1783485           75416   
NJ_0051         5    1654916              1648948           78887   
NJ_0052         5    2807298              2797181           78465   
NJ_0053         5    4384643              4373856           93772   
NJ_0054         5    1844422              1839502           69405   
NJ_0055         5    2190132              2183890          106435   
NJ_0056         5    2799489              2791712          104933   
NJ_0057         5     745198               727419           68637   
NJ_0058         5    3146868              3125059          103023   
NJ_0060         5     978737               976257           56523   
NJ_0061         5    3770196              3760323           75802   
NJ_0062         5     159115               158636           28172   
NJ_0063         5    3581430              3572356          178327   
NJ_0064         5    4964039              4951402           91954   
NJ_0065         5    1796649              1773103           50831   
NJ_0066         5     835816               833685           61751   
NJ_0067         5     993597               974608           77524   
NJ_0068         5    9412347              9388139          187474   
NJ_0069         5   22290334             22236146          248468   
NJ_0070         5   14961759             14926921          301164   
NJ_0071         5     889101               887109           68052   
NJ_0072         5    2341715              2335426           73032   
NJ_0073         5    1433552              1415310           84120   
NJ_0074         5    7854099              7833293          128223   
NJ_0075         5    1009744               761190           56479   
NJ_0076         5     882311               846632           76321   
NJ_0077         5   14442857             14410267          116013   
NJ_0078         5    8557721              8535933           95483   
NJ_0079         5    1399754              1396580          138284   
NJ_0080         5     630119               628747           69109   
NJ_0081         5    2318607              2294005          128507   
NJ_0082         5    8133790              8097476          173446   
NJ_0083         5    2131262              2126143           84477   
NJ_0084         5    2363286              2356677           86809   
NJ_0085         5    2475902              2470464          112703   
NJ_0086         5   11675566             11651298          139996   
NJ_0087         5    7009026              6993284          100563   
NJ_0088         5   12204090             12165858          154044   
NJ_0089         5    3957722              3950461          112799   
NM_3875         5    1113757              1110586           56976   
NM_759          5    2536858              2501458          115954   
NM_761          5     753500               753323           89746   
NM_821          5    1289595              1289300           99947   
NM_839          5     632028               398312           41966   
NM_863          5    1383160              1378601           83393   
NM_863D         5     976271               976084           89355   
NM_864          5     998691               672885           97158   
NM_864D         5     990931               990721          139863   
NM_870          5    1131305              1131074           62898   
NM_871          5    1329276              1325158           68275   
NM_872          5    1041835              1041619           89021   
NM_873          5    1100105              1099872           64867   
NM_877          5     995419               992256           57088   
NM_879          5     853565               853389           56418   
NM_879D         5    2174081              2173587           68398   
NM_885          5    1147719              1144390           57111   
NM_SNBU         5     957172               956963           61917   
NM_SNSA         5    2295401              2294944          234173   
NM_SNSA2        5     311295               311229           63965   

            clusters_hidepth  hetero_est  error_est  reads_consens  
DC_1                     134    0.005308   0.001329            131  
DC_2                   40564    0.016693   0.002049          36789  
DC_3                    8940    0.005010   0.000918           8859  
DC_4                   36749    0.007573   0.001608          35407  
FGXCONTROL              4728    0.001282   0.001397           4673  
ML_1                   80117    0.013519   0.002079          73626  
ML_2                   18778    0.009793   0.001957          16995  
ML_3                    2348    0.009088   0.001667           2300  
ML_4                    6493    0.012995   0.002980           5677  
ML_5                    7615    0.011525   0.002637           6780  
NJ_0003                42230    0.021039   0.003053          35134  
NJ_0005                54567    0.016550   0.001732          47644  
NJ_0006                41121    0.017979   0.001936          35110  
NJ_0007                42620    0.025856   0.002690          34944  
NJ_0008                41098    0.019840   0.002787          34410  
NJ_0009                66563    0.030777   0.003647          52907  
NJ_0010                63185    0.029828   0.003458          50677  
NJ_0011                47935    0.017666   0.001842          41489  
NJ_0012                40850    0.020030   0.002992          34165  
NJ_0015                28088    0.020467   0.002758          23179  
NJ_0026                45043    0.022008   0.002658          37587  
NJ_0029                75645    0.018572   0.001948          66218  
NJ_0030                50793    0.021185   0.002645          42660  
NJ_0031                13496    0.018112   0.002326          11332  
NJ_0032                45314    0.019884   0.001974          38686  
NJ_0033                59051    0.019101   0.001759          51163  
NJ_0034                76757    0.019749   0.002993          64791  
NJ_0035                34729    0.019773   0.002375          28980  
NJ_0036                51253    0.022772   0.002636          42822  
NJ_0037                49383    0.021730   0.002497          41340  
NJ_0038                60158    0.021085   0.001795          51923  
NJ_0039                62243    0.020848   0.002412          52888  
NJ_0040                70184    0.017455   0.001689          60920  
NJ_0041                73180    0.023126   0.002051          61994  
NJ_0045                70981    0.024539   0.002568          59794  
NJ_0047                70956    0.021202   0.002121          60083  
NJ_0048                39186    0.019895   0.001952          33209  
NJ_0049                54130    0.019663   0.002080          46587  
NJ_0051                52604    0.022210   0.002460          44311  
NJ_0052                58132    0.017298   0.001821          50853  
NJ_0053                72883    0.020625   0.002306          62406  
NJ_0054                53049    0.018328   0.001846          46127  
NJ_0055                67244    0.030217   0.003633          53682  
NJ_0056                72063    0.027690   0.002763          59099  
NJ_0057                32991    0.020847   0.002424          27289  
NJ_0058                69527    0.022007   0.002589          58985  
NJ_0060                38635    0.008429   0.002897          35934  
NJ_0061                59500    0.018336   0.001663          51938  
NJ_0062                 5953    0.015460   0.002118           5050  
NJ_0063                86210    0.014303   0.001939          76711  
NJ_0064                67246    0.022482   0.001792          57497  
NJ_0065                18537    0.021576   0.005675          14896  
NJ_0066                32084    0.020127   0.003021          26448  
NJ_0067                30886    0.027449   0.007137          23177  
NJ_0068               104428    0.019692   0.002504          89019  
NJ_0069               108354    0.018701   0.001270          95300  
NJ_0070               136434    0.017976   0.001526         118244  
NJ_0071                37385    0.022200   0.002734          30680  
NJ_0072                55911    0.020707   0.001863          48206  
NJ_0073                50145    0.026550   0.002836          41104  
NJ_0074                88029    0.027581   0.002476          72635  
NJ_0075                29780    0.018548   0.002437          26562  
NJ_0076                38313    0.022905   0.003287          31288  
NJ_0077                77325    0.018676   0.001548          66695  
NJ_0078                71618    0.017639   0.001730          61945  
NJ_0079                60260    0.014759   0.002580          52480  
NJ_0080                24387    0.027919   0.008476          17931  
NJ_0081                72663    0.028462   0.005226          56482  
NJ_0082               118211    0.032335   0.005040          91266  
NJ_0083                60958    0.019438   0.003130          51057  
NJ_0084                64307    0.021354   0.002587          54727  
NJ_0085                62525    0.019904   0.002479          53813  
NJ_0086                87387    0.025240   0.002026          73381  
NJ_0087                72859    0.021598   0.002043          62346  
NJ_0088                93824    0.021661   0.002133          80261  
NJ_0089                45170    0.018329   0.003797          38080  
NM_3875                37545    0.008567   0.001674          35639  
NM_759                 57548    0.012299   0.002059          53875  
NM_761                 37086    0.009551   0.002154          34815  
NM_821                 48760    0.011140   0.002314          45604  
NM_839                 21401    0.009090   0.001479          20921  
NM_863                 51115    0.009782   0.002075          48312  
NM_863D                42406    0.011604   0.002247          39578  
NM_864                 35553    0.008550   0.001873          34091  
NM_864D                49943    0.011054   0.002483          45706  
NM_870                 40324    0.009842   0.002176          37891  
NM_871                 43697    0.009610   0.001973          41202  
NM_872                 46563    0.008602   0.002237          43857  
NM_873                 41092    0.011113   0.001980          38644  
NM_877                 37448    0.009946   0.001902          35367  
NM_879                 34371    0.009148   0.002134          32303  
NM_879D                46517    0.010915   0.002089          43575  
NM_885                 37250    0.008439   0.001619          35422  
NM_SNBU                35248    0.009262   0.002337          33185  
NM_SNSA               117190    0.018197   0.002748         100852  
NM_SNSA2               17634    0.007853   0.002339          16365  


Full stats files
------------------------------------------------
step 1: ./strep_demulti_fastqs/s1_demultiplex_stats.txt
step 2: ./strep_demulti_edits/s2_rawedit_stats.txt
step 3: ./strep_demulti_clust_0.85/s3_cluster_stats.txt
step 4: ./strep_demulti_clust_0.85/s4_joint_estimate.txt
step 5: ./strep_demulti_consens/s5_consens_stats.txt
step 6: None
step 7: Non
```
