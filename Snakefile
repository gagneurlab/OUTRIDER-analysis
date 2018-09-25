configfile: "wbuild.yaml"

htmlOutputPath = config['htmlOutputPath'] if(config['htmlOutputPath'] != None) else 'Output/html'

#
# Global variables combined from config/wbuild.yaml
#
PAPER_METHODS=["pca", "peer", config["FIGURE_IMPLEMENTATION"]]
#, config["FIGURE_ROB_IMPLEMENTATION"]]

include: "./Scripts/Benchmark.snakefile"
include: ".wBuild/wBuild.snakefile"

ruleorder: Scripts_PaperPlots_01d_CreateSimulationNBinomDataSets_R > Scripts_PaperPlots_01c_CreateSimulationNormDataSets_R > Scripts_PaperPlots_03_estimateEncodingDimensions_R
ruleorder: Scripts_PaperPlots_01d_CreateSimulationNBinomDataSets_R > Scripts_PaperPlots_01c_CreateSimulationNormDataSets_R > Scripts_PaperPlots_02_filter_OutriderDataSet_R
    
#
# PAPER FIGURE PIPELINE
#
rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html", "Output/paper_pipeline.done"
    output: touch("Output/all.done")

rule supplement_pdf:
    input: tex='src/latex/supplement_tex/main_supplement.tex', figures=rules.Scripts_PaperPlots_05_FigureX_merge_all_R.input
    output: 'Output/paper_figures/supplement_final.pdf'
    shell: """
        set -x
        cd src/latex/supplement_tex
        pdflatex main_supplement
        bibtex main_supplement
        pdflatex main_supplement
        pdflatex main_supplement
        cp main_supplement.pdf '../../../{output}'
        pwd
        gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/prepress -dNOPAUSE \
            -dQUIET -dBATCH -sOutputFile='../../../{output}.reduced.pdf' '../../../{output}'
    """

rule paper_pipeline:
    input: rules.supplement_pdf.output, rules.Scripts_PaperPlots_06_SummaryStatistics_R.output, 
            rules.Scripts_PaperPlots_05_FigureX_merge_all_R.output
    output: touch("Output/paper_pipeline.done")
# 
# Run full outlier call
# 
rule callOutliers:
    input: expand(config["htmlOutputPath"] + "/datasetSummary/{dataset}.html", dataset=list(set(config["datasets"] + config["GTEx_tissues"])))
    output: touch("Output/run_all_outlier_call_methods.done")

#
# Zscore vs encDim
#
#zscoreImpl=['edPcaRobNfTMax250']
#rule Zscore_vs_EncDim:
#    input: expand("Output/data/zscore_vs_encDim_{dataset}_{impl}.RDS", dataset=["GTEx_not_sun_exposed", "Kremer"], impl=zscoreImpl)
#    output: touch("Output/zscore_vs_encdim.done")
#
#rule Zscore_with_EncDim:
#    input: expand(config["DATADIR"] + "/zscore_with_encDim_{dataset}_{impl}.RDS", dataset=config["datasets"], impl=zscoreImpl)
#    output: touch("Output/zscore_with_encdim.done")


rule runGTEx_variant_enrichment:
    input: expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}_featureset.RDS", dataset=config["GTEx_tissues"])
    output: touch("Output/gtex_enrichment.done")

rule multipleBenchmarkRuns:
    input: "Output/data/fitOutrider/" + config['AE_IMPLEMENTATION'] + "/Kremer_ODS.RDS"
    output: "Output/data/fitOutrider/" + config['AE_IMPLEMENTATION'] + "/Kremer{ID,[0-9]+}_ODS.RDS"
    shell: "ln -s `basename {input}` {output}"
ruleorder: multipleBenchmarkRuns > Scripts_PaperPlots_04_runOutrider_R

rule prepublish:
    shell: """
        rsync -rt {config[htmlOutputPath]} {config[webDir]}
        rsync -rt Output/paper_figures {config[webDir]}/
        """
