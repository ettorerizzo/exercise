import sys, os
import textwrap
import configparser
from configparser import ConfigParser, ExtendedInterpolation
from utils import cmd_init, cmd_out, cmd_out_vcf, _list2input

# Read config File
config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
config.read(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','settings.conf'))
            
def IndelRealigner(core_req, mem_req, in_bam , out_bam, chr, assembly, **tags):
    
    cmd_main = textwrap.dedent(r"""
        outFile={out_bam} && export OUT_FILE=$outFile; 
        
        {config[tools][java]} -Xmx{mem_req}M \
        -jar {config[tools][gatk]} \
        -T RealignerTargetCreator \
        -R {reference} \
        -o {chr}.intervals \
        --known {oneK_indels_vcf} \
        --known {mills_vcf} \
        --num_threads {core_req} \
        -L {chr} \
        {inputs}; 
        
        printf "\n%s RealignerTargetCreator ended.\n" "$(date "+%T %D")" | tee -a /dev/stderr;

        {config[tools][java]} -Xmx{mem_req}M \
        -jar {config[tools][gatk]} \
        -T IndelRealigner \
        -R {reference} \
        -o out.bam \
        -targetIntervals {chr}.intervals \
        -known {oneK_indels_vcf} \
        -known {mills_vcf} \
        -model USE_READS \
        -compress 0 \
        -L {chr} \
        {inputs};
        """.format(config=config,
                   inputs=_list2input(in_bam, '-I '),
                   oneK_indels_vcf=config.get('{assembly}'.format(assembly=assembly),'1kindel_vcf'),
                   mills_vcf=config.get('{assembly}'.format(assembly=assembly),'mills_vcf'),
                   reference=config.get('{assembly}'.format(assembly=assembly),'reference_fasta'),
                   **locals()))

    return cmd_init()+cmd_main+cmd_out()

def HaplotypeCaller(core_req, mem_req, in_bam , out_vcf, chr, assembly, **tags):
    
    cmd_main = textwrap.dedent(r"""
        outFile={out_vcf} && export OUT_FILE=$outFile; 
        
        {config[tools][java]} -Xmx{mem_req}M \
        -jar {config[tools][gatk]} \
        -T HaplotypeCaller \
        -R {reference} \
        -o {out_vcf} \
        -L {chr} \
        -stand_call_conf 30 \
        {inputs};
        """.format(config=config,
                   inputs=_list2input(in_bam, '-I '),
                   dbsnp_vcf=config.get('{assembly}'.format(assembly=assembly),'dbsnp_vcf'),
                   reference=config.get('{assembly}'.format(assembly=assembly),'reference_fasta'),
                   **locals()))

    return cmd_init()+cmd_main+cmd_out_vcf()
