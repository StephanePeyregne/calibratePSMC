# Snakemake Pipeline for Demographic History Reconstruction using PSMC with Heavily Filtered Data (e.g. ancient genomes)

This repository provides instructions for running simulations to correct the effect of filters on the demographic histories estimated with the Pairwise Sequentially Markovian Coalescent (PSMC; Li and Durbin, 2011).
The original workflow was developed by Fabrizio Mafessoni for analyzing high-coverage archaic genomes, as described in Prüfer et al., 2017 and Mafessoni et al., 2020.
The workflow in this repository has been implemented as a snakemake pipeline by Alba Bossoms Mesa, Arev Pelin Sümer, and Stéphane Peyrégne, and was tested with PSMC version 0.6.5-r67.
Example input files are provided in this repository (for chromosomes 21 and 22 from Denisova 3) to run this pipeline.

## 1. Install PSMC

**Clone the Repository**
```
   git clone https://github.com/lh3/psmc.git PSMC
```
**Navigate to the directory**
```
   cd PSMC/
```
**Compile the software**
```
   make
```


## 2. Install scrm (and autoconf)

Before installing scrm, you need to install autoconf.

**Install autoconf**

On Debian/Ubuntu based systems:
```
apt-get install build-essential autoconf autoconf-archive libcppunit-dev
```
On Mac OS:
```
port install automake autoconf autoconf-archive cppunit
```

**Clone the scrm Repository**
```
apt-get install scrm
```

**Build the binary**
```
./bootstrap
make

```


## 3. Prepare input data

1. **Convert VCF Files to PSMC Format**

For each chromosome, convert VCF files to PSMC format and compress the output:
```
   for i in `seq 1 22` ; do python3 scrm2psmc.py -v -g genome_length.txt -m <(zcat filters/chr${i}_mask.bed.gz) <(zcat chr${i}.vcf.gz) | gzip > chr${i}.psmcfa.gz & done
```
The "genome_length.txt" specifies the size of each chromosome, and the "filters/chr${i}_mask.bed.gz" the coordinates of the filters applied to the data that one wishes to replicate.

NOTE: For ancient human genomes, standard filters include a 35-mer mappability filter, GC-dependent coverage filters, repeat filters. As background selection affects demographic inferences, on top of these standard filters, we recommend removing all regions with evidence of more than 10% reduction in genetic variation due to background selection (B-score<900) before running this pipeline.

2. **Concatenate the individual chromosome PSMC files into a single file**
```
   for i in `seq 1 22`; do zcat chr${i}.psmcfa.gz >> Denisova3.psmcfa; done
```
3. **Execute PSMC with default parameters**
```
   PSMC/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o Denisova3.psmc Denisova3.psmcfa
```
4. **Generate simulation code from the PSMC result**
```
   PSMC/utils/psmc2history.pl Denisova3.psmc | PSMC/utils/history2ms.pl -L 2881033286 -u 1.45e-08 -g 29 | sed 's/msHOT-lite/scrm/g' | sed 's/-l//g' > Denisova3_scrm.txt
```
In this command we provide the length of the genome to be simulated (2881033286), the mutation rate (1.45e-08) and the generation time (29)

NOTE: Due to uncertainties in recent colaescent rates, the Ne for the first period is sometimes extremely high, which hampers downstream analyses. This is most likely an artifact and leads to population sizes of 0 in the other periods when converting this history to a scrm command. If you experience similar difficulties, we recommend to manually edit the .hist file to set the problematic population size to that of the next period, and re-running the next steps of the pipeline.

5. **Generate uncalibrated history**
```
   PSMC/utils/psmc2history.pl Denisova3.psmc > Denisova3.history_uncalibrated
```


## 4. Run simulations to calibrate the PSMC results according to the fragmentation of the filtered genome

1. **Edit the Configuration File**

Update config.yaml with your sample name, data paths, and other parameters as shown in the provided example.

2. **Execute the Snakemake pipeline**
```
   snakemake -j 8 -k -R summary_table
```
Change the number of cores indicated in this command (8 cores in this example) as you see fit. PSMC sometimes fails so we use the "-k" --keep-going flag to continue executing the workflow even if some jobs fail, we then re-run the workflow to execute the failed jobs.


## 5. Analyze the results

1. **Determine the correction factors that maximize fit to the real data**
```
   awk 'BEGIN{FS=","; OFS="\t"}{print $1,$2,$3,$4}' calibration/simulations/likelihood_table.csv | sort -k 4 | head -n 1
```
By default, the pipeline explores correction factors ranging from  1.0 to 2.0 in steps of 0.05. Adjust as you see fit. 

Example output where the columns correspond to the name of the sample, the correction factor for theta (scaled mutation rate), the correction factor for rho (scaled recombination rate) and the log-likelihood, in this order.

   Denisova3	1.4	1.6	-1413337.940279

2. **Generate the calibrated demographic history using the best correction factors (e.g., theta = 1.4 and rho = 1.6):**
```
   PSMC/utils/psmc2history.pl calibration/simulations/Denisova3_1.4_1.6.psmc > Denisova3.history_calibrated
```

