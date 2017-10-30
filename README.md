Code in this repo attempts to do open ref BURST using the IMP dataset. Currently about 10% of sequences are not hitting the reference database. With open-ref, we're able to get the number of non-hits down to only 2.2%.

    Reminder of what EMBALMER's output looks like
    Querylabel targetlabel percentid alignmentlength nummismatch numgap startposq endposq startpost endpost evalue bitscore taxonomy 
    T.CS.141_47670	367213	98.909088	275	3	0	1	275	504	779	3	0	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae

1. Grab all sequences that don't hit the database with EMBALMER

    *Teraminx*
    ```bash
        cd /project/flatiron2/pj/imp/all_seq_runs/PROK
        sed 's/\t.*$//' embalmer_output.b6 > seqids.hit.txt
        # transfer to msi
        scp ./seqids.hit.txt vanga015@login.msi.umn.edu:./openref
        cd ../1.trimmed_filtered_fastas
        scp ./combined.seqs.fna vanga015@login.msi.umn.edu:./openref
    ```
    *MSI*
    ```bash
        # run filter.fasta.py on large server (takes too much memory to run interactively)
        qsub -q ram1t filter.fasta.pbs
        # transfer results over back to teraminx
        scp ./combined_seqs_nothit.fna vanga015@teraminx.cs.umn.edu:/project/flatiron2/pj/imp/all_seq_runs/OPENREF
        # note: approx 10% didn't hit
        #    combined_seqs.fna has 20403142 sequences
        #    combined_seqs_nothit.fna has 2130548 sequences
    ```
    
2. On Teraminx, run ninja_filter.exe (use latest version in `/project/flatiron2/sop/bin` - already in path) to deduplicate
    ```bash
    cd $FHOME/imp/all_seq_runs/OPENREF
    mkdir ninja_filter
    /project/flatiron2/sop/bin/ninja_filter_linux ./combined_seqs_nothit.fna D20 D 20
    # "D 2" throws out singletons
    # "D 20" throws out anything that doesn't have at least 20 copies of something
    # add D20 as the prefix for both db file and fa file
    ```
    This will produce two files: D20.db and D20_filt.fa  
    D20_filt.fa = contains all sequences that appeared 20 times or more  
    D20.db = acts like a key to tell you where they appeared  
    Don't need to know this but just in case:  
    * line1: num samples
    * nextlines: list sample names
    * nextlines: for each unique sequence, lists colon-delimited pairs of sample-id-line-number and number of duplicates (92:1 = this sequence appears one time in sample ID found in line (92+2) of D20.db)

    Manually check that this sequence appears correctly in shi7 file
    ```
    grep -B1 '^AATACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGCAGCCGGGTCTGCAAGTCAGATGTGAAATCCATGGGCTCAACCCATGAACTGCATTTGAAACTGTAGATCTTGAGTGTCGGAGGGGCAATCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGATTGCTGGACGATAACTGACGGTGAGGCGCGAAAGTGTGGGGAGCAAACAGGATTAGATACCCGAGTAGTCC' ./combined_seqs_nothit.fna
    ```
    
3. Using sequences that appeared more than 20 times, run all-vs-all embalmer
```
    time emb12 -r D20_filt.fa -q D20_filt.fa  -o /project/flatiron2/pj/imp/all_seq_runs/OPENREF/embalmer_D20.b6 -n -m FORAGE -bs -i 0.935 -f -sa
    # note that to include a more comprehensive list of all-vs-all, decrease -i to .001
```
4. Generate distance matrix from embalmer output
   ```
   Rscript embalmer.to.dm.r -i embalmer_D20.b6 -o D20_dm.txt
   ``` 
5. Calculate optimal number of clusters by generating a bunch of asw scores in parallel (start 4 processes on Teraminx)
```
    Rscript compute.asw.r -i D20_dm.txt -s 2 -e 25 & Rscript compute.asw.r -i D20_dm.txt -s 26 -e 50 & Rscript compute.asw.r -i D20_dm.txt -s 51 -e 75 & Rscript compute.asw.r -i D20_dm.txt -s 76 -e 100
    
    # continue increasing the number of clusters
    Rscript compute.asw.r -i D20_dm.txt -s 101 -e 125 & Rscript compute.asw.r -i D20_dm.txt -s 126 -e 150 & Rscript compute.asw.r -i D20_dm.txt -s 151 -e 175 & Rscript compute.asw.r -i D20_dm.txt -s 176 -e 200
   
   # optimal k = 111!!!!
``` 

6.  Now actually cluster using k=111 on teraminx and save the object to file (to be reopened on local comp). This will output a table of medoid representative IDs.
```
    Rscript ./run.pam.r -i ./D20_dm.txt -k 111 -o pam.111.obj
```
7.  Now let's grab the actual representative sequences and make it our reference db
```
    filter_fasta.py -f D20_filt.fa -s medoid.ids.txt -o rep_seqs.fa
```
#  Steps to only realign failed sequences
8. Now run BURST using these new `rep_seqs.fa` with all of the ones that failed `/project/flatiron2/pj/imp/all_seq_runs/OPENREF/combined_seqs_nothit.fna`
```
            # make your embalmer database files based on the rep_seqs.fa
            burst12 -d -o rep_seqs.edb -a rep_seqs.acc -r rep_seqs.fa -f

            burst12 -r rep_seqs.edb -a rep_seqs.acc -o /dev/null

            burst12 -r rep_seqs.edb -o /dev/null
   
            # Let's try to guess what the taxonomy is for the new DENOVO clusters, and then align
            # Align representative sequences to PROK database at 70%
            burst12 -q rep_seqs.fa -r /project/flatiron2/sop/PROK_170704.edx -a /project/flatiron2/sop/PROK_170704.acx -b /project/flatiron2/sop/PROK_170704.tax -n -m CAPITALIST -bs -i 0.70 -f -o ./embalmer_rep_seqs_70.b6

            # grab the 1st and last columns of the b6 file
            cut -d$'\t' -f 1,13 embalmer_rep_seqs_70.b6 > embalmer_rep_seqs_70_IDs_Taxonomy.txt

            # append this to the PROK taxonomy file
            cat /project/flatiron2/sop/PROK_170704.tax embalmer_rep_seqs_70_IDs_Taxonomy.txt > PROK_IMP_REPSEQS.tax

            # rerun burst with this taxonomy file
            burst12 -r rep_seqs.edx -a rep_seqs.acx -b PROK_IMP_REPSEQS.tax -q /project/flatiron2/pj/imp/all_seq_runs/OPENREF/combined_seqs_nothit.fna -o ./embalmer_output_nothit_tax.b6 -n -m CAPITALIST -bs -i 0.935 -f -sa

            # concatenate this embalmer output with original b6 file, then embalmulate

            cat ../embalmer_output.b6 ./embalmer_output_nothit_tax.b6 > embalmer_output_all_tax.b6

            embalmulate embalmer_output_all_tax.b6 otutable_tax.txt taxatable_tax.txt GGtrim
            # Parsed 19923544 reads [377 samples, 1889 taxa, 2050 refs]. Collating...
            # >>> 19923544/20403142 = 97.6% of all original sequences now hit!!!        
```
9.  Transfer OTU and Taxa files back to local computer, then generate qiime files
Rename otutable_tax.txt to otutable.txt and taxatable_tax.txt to taxatable.txt
Replace all underscores in otutable to spaces
```
    mkdir qiime_files
    bash /Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/embalmer_to_qiime_output.sh otutable_formatted.txt taxatable.txt ./qiime_files /Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/sequences_110716_dimitri_with_PROK_170704/PROK_170704.tre openref
    Rscript /Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/embalmer-post-processing/gabe.burst.summarize.r -i ./qiime_files/taxa0_s2_f.txt
    # manually add "Unknown" for blank taxon name
```        

# Steps to realign all sequences
8. Start by adding the rep_seqs.fa to PROK_170704.fna, create the databases, and use the previously generated taxonomy file (that combined PROK with the 70% aligned rep_set)
```
cat /project/flatiron2/sop/PROK_170704.fna rep_seqs.fa > PROK_170704_REPSEQS.fa

burst12 -d -o PROK_170704_REPSEQS.edb -a PROK_170704_REPSEQS.acc -r PROK_170704_REPSEQS.fa -f

burst12 -r PROK_170704_REPSEQS.edb -a PROK_170704_REPSEQS.acc -o /dev/null

burst12 -r PROK_170704_REPSEQS.edb -o /dev/null
        
# realign ALL sequences (those that hit and didnt hit) against new database + denovo representatives along with the concatenated taxonomy file 
burst12 -r PROK_170704_REPSEQS.edx -a PROK_170704_REPSEQS.acx -b PROK_IMP_REPSEQS.tax -q /project/flatiron2/pj/imp/all_seq_runs/1.trimmed_filtered_fastas/combined_seqs.fna -o ./embalmer_rerun_all_seqs.b6 -n -m CAPITALIST -bs -i 0.935 -f -sa

embalmulate embalmer_rerun_all_seqs.b6 otutable_rerun_all.txt taxatable_rerun_all.txt GGtrim
# Parsed 19964535 reads [377 samples, 1901 taxa, 2044 refs]. Collating...
19964535/20403142 = 97.8% now hit!
```
9. Transfer to local machine. Replace all underscores in otutable to spaces.
```
cd /Users/pvangay/Dropbox/UMN/KnightsLab/Embalmer-Open-Ref/output/with\ all\ seqs\ realigned
bash /Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/embalmer_to_qiime_output.sh otutable_rerun_all.txt taxatable_rerun_all.txt ./ /Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/sequences_110716_dimitri_with_PROK_170704/PROK_170704.tre openref
 
Rscript /Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/embalmer-post-processing/gabe.burst.summarize.r -i ./taxa0_s2_f.txt

```
