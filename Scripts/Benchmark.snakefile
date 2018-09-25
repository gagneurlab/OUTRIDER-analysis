
import os
OUTDIR = os.path.join(config['DATADIR'], config['BENCHMARK_DIR'])
RECALL_FILES = expand(OUTDIR + 
        "/nsamples={nsamples}/inj={inj}/{inj_value}/correction={correction}" +
        "/q={q}/pmethod={pmethod}/{dataset}_plot.tsv",
    nsamples=config['N_samples'],
    inj=config['inj'],
    inj_value=config['inj_value'],
    correction=config['correction'],
    q=config['Qs'],
    pmethod=config['FDR_METHOD'],
    dataset=config['datasets'])

#
# GTEx variant enrichment
#
rule extractAllRareAnnoFromGTEx: 
    input:
        vcfanno = config["GTEXVCF_ANNO"]
    output:
        subvcfannotbi = config["DATADIR"] + "/GTEx_variant_enrichment/rareAllVariants.vcf.gz.tbi"
    shell:
        """
            vcfAnnoFileTbi={output.subvcfannotbi}
            maxAF=0.05
            
            vcfAnnoFile=${{vcfAnnoFileTbi%.tbi}}
            bcftools view --max-af $maxAF --max-alleles 2 {input.vcfanno} | \
                grep -P '(^#)|(\|(protein_coding)|(lincRNA)\|)' | \
                bgzip -c > $vcfAnnoFile
            tabix $vcfAnnoFile
        """
        
rule extractModeratNHighAnnoFromGTEx: 
    input:
        vcfanno = config["GTEXVCF_ANNO"]
    output:
        subvcfannotbi = config["DATADIR"] + "/GTEx_variant_enrichment/rareModerateNHighVariants.vcf.gz.tbi"
    shell:
        """
            vcfAnnoFileTbi={output.subvcfannotbi}
            maxAF=0.05
            
            vcfAnnoFile=${{vcfAnnoFileTbi%.tbi}}
            bcftools view --max-af $maxAF --max-alleles 2 {input.vcfanno} | \
                grep -P '(^#)|(\|(MODERATE|HIGH)\|)' | \
                bgzip -c > $vcfAnnoFile
            tabix $vcfAnnoFile
        """


rule extractGenoTypeFromGTEx:
    input:
        subvcfannotbi = config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}.vcf.gz.tbi",
        vcfgt =   config["GTEXVCF_GT"]
    output:
        subvcfgttbi =  config["DATADIR"]  + "/GTEx_variant_enrichment/{dataset}_gt.vcf.gz.tbi"
    shell:
        """
            set -x
            vcfAnnoFileTbi={input.subvcfannotbi}
            vcfGTFileTbi={output.subvcfgttbi}
            maxAF=0.05
            
            vcfAnnoFile=${{vcfAnnoFileTbi%.tbi}}
            vcfGTFile=${{vcfGTFileTbi%.tbi}}
            
            bcftools view --regions-file $vcfAnnoFile --max-af $maxAF {input.vcfgt} | \
                bcftools sort --max-mem 10G --temp-dir '{input.subvcfannotbi}.tmp/' | bgzip -c > $vcfGTFile
            tabix $vcfGTFile
        """


#
# RECALL pipeline
#
rule Recall:
    input: RECALL_FILES, 
            expand(config["htmlOutputPath"] + 
                    "/RecallBenchmark/da_{dataset}_co_{correction}_ggplots.RDS", 
                    correction=config['correction'], dataset=config['datasets']),
            expand(config["htmlOutputPath"] + 
                    "/RecallBenchmark/da_{dataset}_all_ggplots.RDS",
                    dataset=config['datasets'])
    


rule injection:
    """Generate the injection mask"""
    input:
        counts = config['DATADIR'] + '/fitOutrider/' + config['AE_IMPLEMENTATION'] + '/{dataset}_ODS.RDS',
        scriptFile = "src/r/AnalysisPipeline/injectionMask.R"
    output:
        true_outliers = OUTDIR + "/nsamples={nsamples}/{dataset}_outlier_mask.RDS",
    script:
        "../src/r/AnalysisPipeline/injectionMask.R"

rule preprocess:
    """Generate the injected count table"""
    input:
        true_outliers = OUTDIR + "/nsamples={nsamples}/{dataset}_outlier_mask.RDS",
        scriptFile = "src/r/AnalysisPipeline/injectOutliers.R"
    output:
        inj_counts = OUTDIR + "/nsamples={nsamples}/inj={inj}/{inj_value}/{dataset}_injectedcounts.RDS"
    script:
        "../src/r/AnalysisPipeline/injectOutliers.R"

rule fit_outrider:
    """Runs the outrider methods"""
    input:
        inj_counts = OUTDIR + "/nsamples={nsamples}/inj={inj}/{inj_value}/{dataset}_injectedcounts.RDS",
        bestQ = config['DATADIR'] + '/fitOutrider/' + config['AE_IMPLEMENTATION'] + '/{dataset}_ODS.RDS',
        KremerCovars = 'Data/rawdata/Kremer_sample_annotation.tsv',
        GTExCovars = 'Data/rawdata/V6P/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt',
        scriptFile = "src/r/AnalysisPipeline/fitOUTRIDER.R"
    output:
        outrider_fit = OUTDIR + "/nsamples={nsamples}/inj={inj}/{inj_value}/correction={correction}/q={q}/{dataset}_OUTRIDERfit.RDS"
    threads: 15
    script:
        "../src/r/AnalysisPipeline/fitOUTRIDER.R"
    

rule run_outrider_pVal:
    """Runs the outrider pVal computation"""
    input:
        outrider_fit = OUTDIR + "/nsamples={nsamples}/inj={inj}/{inj_value}/correction={correction}/q={q}/{dataset}_OUTRIDERfit.RDS",
        scriptFile = "src/r/AnalysisPipeline/computePvalOUTRIDER.R"
    output:
        outrider_pval = OUTDIR + "/nsamples={nsamples}/inj={inj}/{inj_value}/correction={correction}/q={q}/pmethod={pmethod}/{dataset}_OUTRIDERpVal.RDS"
    threads: 4
    script:
        "../src/r/AnalysisPipeline/computePvalOUTRIDER.R"


rule evalOUTRIDER:
    """ranks counts"""
    input:
        outrider_pval = OUTDIR + "/nsamples={nsamples}/inj={inj}/{inj_value}/correction={correction}/q={q}/pmethod={pmethod}/{dataset}_OUTRIDERpVal.RDS",
        scriptFile = "src/r/AnalysisPipeline/evalOUTRIDER.R"
    output:
        plot_table = OUTDIR + "/nsamples={nsamples}/inj={inj}/{inj_value}/correction={correction}/q={q}/pmethod={pmethod}/{dataset}_plot.tsv"
    script:
        "../src/r/AnalysisPipeline/evalOUTRIDER.R"


