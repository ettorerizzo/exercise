import sys, os
import textwrap
import configparser
from configparser import ConfigParser, ExtendedInterpolation
from utils import cmd_init, cmd_out, _list2input, awk_cmd

# Read config File
config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
config.read(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','settings.conf'))

def Fastqc(core_req,mem_req,in_txts,out_dir,out_subdir,out_summary):
    
    cmd_main = textwrap.dedent(r"""
        outDIR={out_dir} && export OUT_DIR=$outDIR; 
           
        {config[tools][fastqc]} \
        -d $tmpDir \
        -f fastq \
        --nogroup \
        --casava \
        -t 4 \
        -o $tmpDir \
        {in_txts[0]} \
        {in_txts[1]} \
        && \
        unzip $tmpDir/\*.zip -d $tmpDir \
        && \
        rm $tmpDir/*.zip \
        && \
        cp -R $tmpDir/* $outDIR \
        && \
        cp -R $tmpDir/* {out_subdir} \
        && \
        find {out_subdir} -type f -name 'summary.txt' | xargs cat > {out_summary} \
        && \
        for file in $(find {out_subdir} -type f -name 'fastqc_data.txt');
        do {awk} $file >> {out_summary};
        done;
        """.format(config=config,
                   awk=awk_cmd(),
                   **locals()))

    return cmd_init()+cmd_main

def FastqcStatistic(core_req, mem_req, in_summary, in_dir, out_dir, out_html):
    
    cmd_main = textwrap.dedent(r"""
        outDIR={out_dir} && export OUT_DIR=$outDIR; 
        
        cp -R {in_dir[0]}* $outDIR \
        && \
        cat \
        {inputs} > $outDIR/summary.txt \
        && \
        sort -k 2 $outDIR/summary.txt > $outDIR/summary.sort.txt \
        && \
        java -jar {config[tools][fastqc_html]} -i $outDIR/summary.sort.txt -o {out_html}
        """.format(config=config,
                   inputs=_list2input(in_summary,""),
                   **locals()))

    return cmd_init()+cmd_main