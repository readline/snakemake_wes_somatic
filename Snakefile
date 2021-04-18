from os.path import join
import os
import pandas as pd
from scripts.Load import somsamplesheet

snakedir = os.getcwd()
print(snakedir)
configfile: 'config.yaml'
print(config)
somdic, ssdic = somsamplesheet(config['samplesheet'])
workdir: config['workdir']
#sbcmd="sbatch --cpus-per-task={threads} --mem={cluster.mem}"
#sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
#sbcmd+=" --out={cluster.out} {cluster.extra}"
print(somdic)

rule all:
    input:
        lambda wildcards: ["12.SV.manta/{}/results/variants/candidateSmallIndels.vcf.gz".format(som) for som in somdic],
        lambda wildcards: ["11.Somatic.Strelka/{}/results/variants/somatic.snvs.vcf.gz".format(som) for som in somdic],
        lambda wildcards: ["11.Somatic.Mutect/{}/{}.mutect2.flt.vcf".format(som, som) for som in somdic],
        lambda wildcards: ["11.Somatic.Lofreq/{}/{}.lofreqsomatic_final.snvs.vcf.gz".format(som, som) for som in somdic],
        lambda wildcards: ["11.Somatic.MuSE/{}/{}.MuSE.vcf.gz".format(som, som) for som in somdic],
        lambda wildcards: ["11.Somatic.Varscan2/{}/{}.snp.LOH.vcf".format(som, som) for som in somdic],

        lambda wildcards: ["11.Somatic.MutectSS/{}/{}.mutect2.flt.vcf".format(som, som) for som in ssdic],
        
rule som_manta:
    input:
        bam1=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][0], somdic[wildcards.som][0]),
        bam0=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][1], somdic[wildcards.som][1]),
    output:
        folder=directory("12.SV.manta/{som}"),
        result="12.SV.manta/{som}/results/variants/candidateSmallIndels.vcf.gz",
    log:
        out = snakedir+"/logs/D01.som_manta/{som}.o",
        err = snakedir+"/logs/D01.som_manta/{som}.e",
    threads:  16
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:10 ',
    shell:
        '''module load {config[modules][manta]}
        configManta.py \
          --normalBam {input.bam0} \
          --tumorBam {input.bam1} \
          --referenceFasta {config[references][fasta]} \
          --runDir {output.folder} \
          --callRegions {config[references][flankbedgz]} >{log.out} 2>{log.err}
        {output.folder}/runWorkflow.py -m local -j 16  >>{log.out} 2>>{log.err} '''

rule som_strelka:
    input:
        bam1=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][0], somdic[wildcards.som][0]),
        bam0=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][1], somdic[wildcards.som][1]),
        manta="12.SV.manta/{som}/results/variants/candidateSmallIndels.vcf.gz",
    output:
        folder=directory("11.Somatic.Strelka/{som}"),
        result="11.Somatic.Strelka/{som}/results/variants/somatic.snvs.vcf.gz",
    log:
        out = snakedir+"/logs/D02.som_strelka/{som}.o",
        err = snakedir+"/logs/D02.som_strelka/{som}.e",
    threads:  16
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:10 ',
    shell:
        '''module load {config[modules][strelka]}
        configureStrelkaSomaticWorkflow.py \
          --tumorBam {input.bam1} \
          --normalBam {input.bam0} \
          --referenceFasta {config[references][fasta]} \
          --runDir {output.folder} \
          --callRegions {config[references][flankbedgz]} \
          --indelCandidates {input.manta} \
          --outputCallableRegions >{log.out} 2>{log.err}
        {output.folder}/runWorkflow.py -m local -j 16  >>{log.out} 2>>{log.err}'''
        
        
rule som_mutect:
    input:
        bam1=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][0], somdic[wildcards.som][0]),
        bam0=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][1], somdic[wildcards.som][1]),
    params:
        normal=lambda wildcards: somdic[wildcards.som][1],
    output:
        vcf1="11.Somatic.Mutect/{som}/{som}.mutect2.vcf",
        f1r2=temp("11.Somatic.Mutect/{som}/{som}.f1r2.tar.gz"),
        model=temp("11.Somatic.Mutect/{som}/{som}.read-orientation-model.tar.gz"),
        pile=temp("11.Somatic.Mutect/{som}/{som}.getpileupsummaries.table"),
        seg=temp("11.Somatic.Mutect/{som}/{som}.segments.table"),
        cont=temp("11.Somatic.Mutect/{som}/{som}.calculatecontamination.table"),
        vcf2="11.Somatic.Mutect/{som}/{som}.mutect2.flt.vcf",
    log:
        out = snakedir+"/logs/D03.som_mutect/{som}.o",
        err = snakedir+"/logs/D03.som_mutect/{som}.e",
    threads:  8
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:10 ',
    shell:
        '''module load {config[modules][gatk]}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          Mutect2 \
          -R {config[references][fasta]} \
          -L {config[references][flankbed]} \
          -I {input.bam1} \
          -I {input.bam0} \
          -normal {params.normal} \
          --panel-of-normals {config[references][mutectpon]}\
          -O {output.vcf1} \
          -germline-resource {config[references][afonlygnomad]} \
          --f1r2-tar-gz {output.f1r2} >{log.out} 2>{log.err}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          LearnReadOrientationModel \
          -O {output.model} \
          -I {output.f1r2} >>{log.out} 2>>{log.err}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          GetPileupSummaries \
          -I {input.bam1} \
          -V {config[references][exaccommon]} \
          -L {config[references][flankbed]} \
          -O {output.pile} >>{log.out} 2>>{log.err}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          CalculateContamination \
          -I {output.pile} \
          -tumor-segmentation {output.seg} \
          -O {output.cont} >>{log.out} 2>>{log.err}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          FilterMutectCalls \
          --reference {config[references][fasta]} \
          -V {output.vcf1} \
          --tumor-segmentation {output.seg} \
          --ob-priors {output.model} \
          -O {output.vcf2} >>{log.out} 2>>{log.err} '''

rule som_mutectss:
    input:
        bam1=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(ssdic[wildcards.som][0], ssdic[wildcards.som][0]),
    params:
        normal=lambda wildcards: ssdic[wildcards.som][1],
    output:
        vcf1="11.Somatic.MutectSS/{som}/{som}.mutect2.vcf",
        f1r2=temp("11.Somatic.MutectSS/{som}/{som}.f1r2.tar.gz"),
        model=temp("11.Somatic.MutectSS/{som}/{som}.read-orientation-model.tar.gz"),
        pile=temp("11.Somatic.MutectSS/{som}/{som}.getpileupsummaries.table"),
        seg=temp("11.Somatic.MutectSS/{som}/{som}.segments.table"),
        cont=temp("11.Somatic.MutectSS/{som}/{som}.calculatecontamination.table"),
        vcf2="11.Somatic.MutectSS/{som}/{som}.mutect2.flt.vcf",
    log:
        out = snakedir+"/logs/D07.som_mutectss/{som}.o",
        err = snakedir+"/logs/D07.som_mutectss/{som}.e",
    threads:  8
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:10 ',
    shell:
        '''module load {config[modules][gatk]}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          Mutect2 \
          -R {config[references][fasta]} \
          -L {config[references][flankbed]} \
          -I {input.bam1} \
          --panel-of-normals {config[references][mutectpon]}\
          -O {output.vcf1} \
          -germline-resource {config[references][afonlygnomad]} \
          --f1r2-tar-gz {output.f1r2} >{log.out} 2>{log.err}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          LearnReadOrientationModel \
          -O {output.model} \
          -I {output.f1r2} >>{log.out} 2>>{log.err}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          GetPileupSummaries \
          -I {input.bam1} \
          -V {config[references][exaccommon]} \
          -L {config[references][flankbed]} \
          -O {output.pile} >>{log.out} 2>>{log.err}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          CalculateContamination \
          -I {output.pile} \
          -tumor-segmentation {output.seg} \
          -O {output.cont} >>{log.out} 2>>{log.err}
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          FilterMutectCalls \
          --reference {config[references][fasta]} \
          -V {output.vcf1} \
          --tumor-segmentation {output.seg} \
          --ob-priors {output.model} \
          -O {output.vcf2} >>{log.out} 2>>{log.err} '''
        
rule som_lofreq:
    input:
        bam1=lambda wildcards: "02.Alignment/Lofreq/{}/{}.li.bam".format(somdic[wildcards.som][0], somdic[wildcards.som][0]),
        bam0=lambda wildcards: "02.Alignment/Lofreq/{}/{}.li.bam".format(somdic[wildcards.som][1], somdic[wildcards.som][1]),
    params:
        prefix="11.Somatic.Lofreq/{som}/{som}.lofreq",
    output:
        vcf1="11.Somatic.Lofreq/{som}/{som}.lofreqsomatic_final.snvs.vcf.gz",
    log:
        out = snakedir+"/logs/D04.som_lofreq/{som}.o",
        err = snakedir+"/logs/D04.som_lofreq/{som}.e",
    threads:  16
    resources:
        mem  = '64g',
        extra = ' --gres=lscratch:20 ',
    shell:
        '''module load {config[modules][lofreq]}
        lofreq somatic \
          -n {input.bam0} \
          -t {input.bam1} \
          -o {params.prefix} \
          -f {config[references][fasta]} \
          --threads 16 \
          --call-indels \
          -d {config[references][snp138]} >{log.out} 2>{log.err} '''

rule som_muse:
    input:
        bam1=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][0], somdic[wildcards.som][0]),
        bam0=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][1], somdic[wildcards.som][1]),
    params:
        prefix="11.Somatic.MuSE/{som}/{som}",
    output:
        vcf="11.Somatic.MuSE/{som}/{som}.MuSE.vcf.gz",
    log:
        out = snakedir+"/logs/D05.som_muse/{som}.o",
        err = snakedir+"/logs/D05.som_muse/{som}.e",
    threads:  4
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:20 ',
    shell:
        '''module load {config[modules][muse]} {config[modules][samtools]}
        MuSE call \
          -f {config[references][fasta]} \
          {input.bam1} \
          {input.bam0} \
          -l {config[references][flankbed]} \
          -O {params.prefix} > {log.out} 2> {log.err}
        MuSE sump \
          -I {params.prefix}.MuSE.txt \
          -G \
          -D {config[references][snp138]} \
          -O {params.prefix}.MuSE.vcf >> {log.out} 2>> {log.err}
        bgzip {params.prefix}.MuSE.vcf
        tabix -p vcf {output.vcf}'''
        

rule som_mpileup:
    input:
        bam1=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][0], somdic[wildcards.som][0]),
        bam0=lambda wildcards: "02.Alignment/Level3/{}/{}.BQSR.bam".format(somdic[wildcards.som][1], somdic[wildcards.som][1]),
    output:
        pile=temp("11.Somatic.Varscan2/{som}/{som}.mpileup.gz"),
        pile0="11.Somatic.Varscan2/{som}/{som}.0.mpileup",
        pile1="11.Somatic.Varscan2/{som}/{som}.1.mpileup",
    log:
        out = snakedir+"/logs/D06.som_varscan/{som}.pileup.o",
        err = snakedir+"/logs/D06.som_varscan/{som}.pileup.e",
    threads:  4
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:20 ',
    shell:
        '''module load {config[modules][samtools]}
        samtools mpileup \
          -q 1 \
          -f {config[references][fasta]} \
          -l {config[references][flankbed]} \
          {input.bam0} \
          {input.bam1} 2>{log.err} | gzip -c - > {output.pile}
        zcat {output.pile}|cut -f 1,2,3,4,5,6 > {output.pile0}
        zcat {output.pile}|cut -f 1,2,3,7,8,9 > {output.pile1}'''
        
rule som_varscan:
    input:
        pile0="11.Somatic.Varscan2/{som}/{som}.0.mpileup",
        pile1="11.Somatic.Varscan2/{som}/{som}.1.mpileup",
    params:
        prefix="11.Somatic.Varscan2/{som}/{som}",
    output:
        loh="11.Somatic.Varscan2/{som}/{som}.snp.LOH.vcf",
    log:
        out = snakedir+"/logs/D06.som_varscan/{som}.var.o",
        err = snakedir+"/logs/D06.som_varscan/{som}.var.e",
    threads:  4
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:20 ',
    shell:
        '''module load {config[modules][varscan]} {config[modules][samtools]}
        varscan \
          somatic \
          {input.pile0} \
          {input.pile1} \
          {params.prefix} \
          --min-coverage 8 \
          --min-coverage-normal 8 \
          --min-coverage-tumor 6 \
          --min-var-freq 0.10 \
          --min-freq-for-hom 0.75 \
          --output-vcf 1 > {log.out} 2>> {log.err}
        varscan \
          processSomatic \
          {params.prefix}.snp.vcf \
          --min-tumor-freq 0.1 \
          --max-normal-freq 0.05 \
          --p-value 0.07 >> {log.out} 2>> {log.err}
        varscan \
          processSomatic \
          {params.prefix}.indel.vcf \
          --min-tumor-freq 0.1 \
          --max-normal-freq 0.05 \
          --p-value 0.07 >> {log.out} 2>> {log.err} &&\
        rm {input.pile0} {input.pile1}'''

        



