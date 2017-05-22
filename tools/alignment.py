import sys, os
import textwrap
import configparser
from configparser import ConfigParser, ExtendedInterpolation
from utils import cmd_init, cmd_out, _list2input

# Read config File
config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
config.read(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','settings.conf'))

def Alignment(core_req,mem_req,in_txts,out_bam,chunk,sample_name,rgid,platform,platform_unit,procedure,library,assembly,result_dir):
    cmd_main = textwrap.dedent(r"""
        outFile={out_bam} && export OUT_FILE=$outFile; 
         
        {config[tools][bwa]} mem \
        -M \
        -t {core_req} \
        -R "@RG\tID:{rgid}\tLB:{library}\tSM:{sample_name}\tPL:{platform}\tPU:{platform_unit}" \
        {reference} \
        {in_txts[0]} \
        {in_txts[1]} \
        | \
        {config[tools][java]} -Xmx{mem_req}M -jar {config[tools][picard_dir]}/AddOrReplaceReadGroups.jar \
        INPUT=/dev/stdin \
        OUTPUT=/dev/stdout \
        RGID={rgid} \
        RGLB={library} \
        RGSM={sample_name} \
        RGPL={platform} \
        RGPU={platform_unit} \
        COMPRESSION_LEVEL=0 \
        | \
        {config[tools][java]} -Xmx{mem_req}M -jar {config[tools][picard_dir]}/CleanSam.jar \
        INPUT=/dev/stdin \
        OUTPUT=/dev/stdout \
        VALIDATION_STRINGENCY=SILENT \
        COMPRESSION_LEVEL=0 \
        | \
        {config[tools][java]} -Xmx{mem_req}M -jar {config[tools][picard_dir]}/SortSam.jar \
        INPUT=/dev/stdin \
        OUTPUT=out.bam \
        SORT_ORDER=coordinate \
        CREATE_INDEX=True
        """.format(config=config,
                   reference=config.get('{assembly}'.format(assembly=assembly),'reference_fasta'),
                   **locals()))

    return cmd_init()+cmd_main+cmd_out()