nextflow.enable.dsl=2
//params.input = "/home/centos/nf-inmputation/testData/chr*.vcf.gz"
//vcf_ch = Channel.fromPath(params.input)
params.input = "/chr*.vcf.gz"
params.vcffiles = params.vcfdir + params.input
vcf_ch = Channel.fromPath(params.vcffiles)
params.fasta_ref = '/home/centos/nf-imputation/testData/reference/Homo_sapiens_assembly38.fasta'

process tabix {
  publishDir "./"
  input:
  path vcf

  output:
  path "${vcf}.tbi"

  script:
  """
  #!/bin/bash
  echo $PATH
  sleep 300 
  module load bcftools
  tabix -f ${vcf}
  """
}

process qc1 {
  publishDir "./"
  input:
  path vcf

  output:
  path "${vcf.SimpleName}_qc1.bcf"

  script:
  """
  #!/bin/bash
  module load bcftools
  bcftools norm -m- -Oz -o ${vcf.SimpleName}_qc1.bcf ${vcf}
  """
}


process qc2 {
  publishDir "./"
  input:
  path bcf

  output:
  path "${bcf.SimpleName}_qc2.bcf"

  script:
  """
  module load bcftools
  bcftools query -f '%CHROM\\t%POS\\t%REF%ALT\\n' ${bcf} > ${bcf.SimpleName}.pos
  awk '{if ( \$3!="AC" && \$3!="CA" && \$3!="CT" && \$3!="TC" && \$3!="TG" && \$3!="GT" && \$3!="AG" && \$3!="GA" ) print \$1"\t"\$2}' ${bcf.SimpleName}.pos > ${bcf.SimpleName}.atgc
  bcftools view -T ^${bcf.SimpleName}.atgc -Ob -o ${bcf.SimpleName}_qc2.bcf ${bcf}
  """
}

process qc3 {
  publishDir "./"
  input:
  path bcf

  output:
  path "${bcf.SimpleName}_qc3.bcf"

  script:
  """
  module load bcftools
  bcftools norm -m+ ${bcf} | bcftools view -v snps -m2 -M2 | bcftools annotate -x FORMAT,INFO | bcftools +fill-tags -- -t all | bcftools view -i 'AC>0' | bcftools view -i 'F_MISSING<0.05' | python3 -c "import sys;[print (line.strip().replace('|','/').replace('1/0','0/1')) if line.startswith('#')==False else print (line.strip()) for line in sys.stdin]" | bcftools view -Ob -o ${bcf.SimpleName}_qc3.bcf && tabix -f ${bcf.SimpleName}_qc3.bcf
  """
}

process qc4 {
  publishDir "./"
  input:
  path bcf

  output:
  path "${bcf.SimpleName}_qc4.bcf"

  script:
  """
  module load bcftools any/plink2
  plink2 --bcf ${bcf} --missing --out ${bcf}
  cat ${bcf}.smiss | sed 1d | awk '{if(\$4 > 0.05) print \$1}' > ${bcf}_mind.txt
  bcftools view -S ^${bcf}_mind.txt ${bcf} | bcftools annotate -x FORMAT | bcftools +fill-tags -- -t all | bcftools view -Ob -o ${bcf.SimpleName}_qc4.bcf && tabix -f ${bcf.SimpleName}_qc4.bcf
  """
}

process qc5 {
  publishDir "./", mode: 'copy', overwrite: true
  input:
  path bcf

  output:
  path "${bcf.SimpleName}_qc5.bcf"
  path "${bcf.SimpleName}_qc5.bcf.csi"
  path "${bcf.SimpleName}_qc5.log"

  script:
  """
  module load bcftools 
  bcftools +fixref ${bcf} -Ob -o ${bcf.SimpleName}_qc5.bcf -- -f ${params.fasta_ref} -m flip > ${bcf.SimpleName}_qc5.log 2>&1
  tabix -f ${bcf.SimpleName}_qc5.bcf
  """
}

workflow {
  vcf_ch.view()
  tabix (vcf_ch)
  qc1 (vcf_ch)
  qc2 (qc1.out)
  qc3 (qc2.out)
  qc4 (qc3.out)
  qc5 (qc4.out)
}

