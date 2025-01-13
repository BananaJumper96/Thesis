######################################################################################################
# Workflow for bin-level metabolic modeling across multiple samples.
#
# Usage:
#   1. Update your config.yaml.
#   2. Ensure each sample directory under SAMPLES_PATH has the same structure:
#         sampleX/
#           bins.contigs.gz
#           annotation.filt.contig2ID.tsv.gz
#           annotation.faa.gz
#           bins.stats.gz (optional)
#   3. Run snakemake as usual with --use-conda, etc.
#
# This Snakefile includes:
#   - A checkpoint MakeBinProteinFastas that extracts per-bin protein FASTAs for each sample.
#   - Rule carveme to build GEMs for each bin.
#   - Rules ECvis, smetana, memote to operate once per sample,
#       aggregating all bin-level GEMs in that sample.
#   - Rule modelVis to get overview over all samples.
######################################################################################################

import os
import glob

configfile: "config.yaml"

OUT_PATH = config["out_path"]
SAMPLES_PATH = config["samples_path"]

ROOTDIR = workflow.basedir
SCRIPTS_DIR = os.path.join(ROOTDIR, "Scripts")
ENV_DIR = os.path.join(ROOTDIR, "Envs")
DB_DIR = os.path.join(ROOTDIR, "Databases")


# Ensure SAMPLES_PATH ends with a slash to match directories
if not SAMPLES_PATH.endswith('/'):
    SAMPLES_PATH += '/'

# Use glob to find all directories in SAMPLES_PATH
SAMPLES = [os.path.basename(d) for d in glob.glob(f"{SAMPLES_PATH}*") if os.path.isdir(d)]

######################################################################################################
# Helpers: Collect all bin-level GEMs for a given sample or for all samples
######################################################################################################
def gatherbinmodels(wildcards):
    """
    Returns all .xml GEMs for the given sample, after per-bin FASTAs are generated.
    This ensures each rule that wants a single set of GEMs per sample gets them all at once.
    """
    # Wait for the checkpoint output
    ckpt_out = checkpoints.MakeBinProteinFastas.get(sample=wildcards.sample).output.protbinfastaDir
    
    # Extract bin IDs
    bin_pattern = os.path.join(ckpt_out, "{binID}.faa")
    binIDs = glob_wildcards(bin_pattern).binID
    
    # Expand to return paths for all GEMs
    return expand(
        os.path.join(OUT_PATH, "{sample}", "GEMs", "{binID}.xml"),
        sample=wildcards.sample,
        binID=binIDs
    )

def gather_all_models(wildcards):
    """
    Gathers ALL .xml GEMs across ALL samples using the MakeBinProteinFastas checkpoint
    """
    all_xml_paths = []
    for sample in SAMPLES:
        # Wait for the checkpoint output per-sample
        ckpt_out = checkpoints.MakeBinProteinFastas.get(sample=sample).output.protbinfastaDir
        
        # Extract bin IDs
        bin_pattern = os.path.join(ckpt_out, "{binID}.faa")
        binIDs = glob_wildcards(bin_pattern).binID
        
        # Expand to return paths for all GEMs in that sample
        xmls_this_sample = expand(
            os.path.join(OUT_PATH, "{sample}", "GEMs", "{binID}.xml"),
            sample=sample,
            binID=binIDs
        )
        
        all_xml_paths.extend(xmls_this_sample)

    return all_xml_paths
def gather_all_ecfiles(wildcards):
    """
    Gathers ALL .ec files across ALL samples using the ECvis rule
    """
    all_ec_paths = []
    for sample in SAMPLES:
        # Wait for the ECvis output per-sample
        ec_dir = rules.ECvis.get(sample=sample).output.ec_dir
        
        # Extract .ec files
        ec_pattern = os.path.join(ec_dir, "*.ec")
        ec_files = glob.glob(ec_pattern)
        
        all_ec_paths.extend(ec_files)

    return all_ec_paths

######################################################################################################
# Master rule: require final outputs for all samples
######################################################################################################
rule all:
    input:
        # modelVis PDF + GEMs.stats
        # os.path.join(OUT_PATH, "Stats", "modelVis.pdf"),
        # os.path.join(OUT_PATH, "GEMs", "GEMs.stats"),
        # SMETANA single TSV per sample
        # expand(os.path.join(OUT_PATH, "{sample}", "SMETANA", "{sample}_detailed.tsv"), sample=SAMPLES),
        # EC directory per sample
        # expand(os.path.join(OUT_PATH, "{sample}", "ecfiles"), sample=SAMPLES),
        # Memote summary per sample
        # expand(os.path.join(OUT_PATH, "{sample}", "memote", "{sample}.html"), sample=SAMPLES),
        tsv = directory(os.path.join(OUT_PATH, "Stats/GEMs_vs_BinCompleteness"))
    message:
        """
        Final per-sample outputs for bin-level metabolic modeling.
        """

######################################################################################################
# Checkpoint: create per-bin protein FASTAs for each sample
######################################################################################################
checkpoint MakeBinProteinFastas:
    input:
        cont2bin      = os.path.join(SAMPLES_PATH, "{sample}", "bins.contigs.gz"),
        cont2orf      = os.path.join(SAMPLES_PATH, "{sample}", "annotation.filt.contig2ID.tsv.gz"),
        protein_fasta = os.path.join(SAMPLES_PATH, "{sample}", "annotation.faa.gz")
    output:
        protbinfastaDir = directory(os.path.join(OUT_PATH, "{sample}", "protein_bins"))
    params:
        bin_quali = lambda wildcards: os.path.join(SAMPLES_PATH, wildcards.sample, f"{config["bin_quali_file_name"]}"),
        min_comp = config["min_comp"],
        max_cont = config["max_cont"]
    log:
        os.path.join(OUT_PATH, "{sample}", "logs", "MakeBinProteinFastas.log")
    conda:
        os.path.join(ENV_DIR, "metaGEMmod.yaml")
    shell:
        """
        mkdir -p $(dirname {log})
        python {SCRIPTS_DIR}/getbinprotfastas.py \
            --cont2bin {input.cont2bin} \
            --cont2orf {input.cont2orf} \
            --protein_fasta {input.protein_fasta} \
            --output_dir {output.protbinfastaDir} \
            --bin_stats {params.bin_quali} \
            --max_cont {params.max_cont} \
            --min_comp {params.min_comp} &> {log}
        """

######################################################################################################
# Rule: generate GEMs via CarveMe (one .xml per bin)
######################################################################################################

rule carveme:
    input:
        bin   = os.path.join(OUT_PATH, "{sample}", "protein_bins", "{binID}.faa"),
        media = os.path.join(DB_DIR, "media_db.tsv")
    output:
        os.path.join(OUT_PATH, "{sample}", "GEMs", "{binID}.xml")
    params:
        binID = lambda w: w.binID,
        carveme_medium  = config["carvemeMedium"],
    conda:
        os.path.join(ENV_DIR, "carveme.yaml")
    log:
        os.path.join(OUT_PATH, "{sample}", "logs", "carveme", "carveme_{binID}.log")
    shell:
        """
        # Create necessary directories for output and logs
        mkdir -p "$(dirname "{output}")"
        mkdir -p "$(dirname "{log}")"

        # Create a unique temporary directory in /dev/shm
        tmpdir=$(mktemp -d -p /dev/shm carveme_XXXXXX)
        echo "Temporary working directory: $tmpdir" >> "{log}"

        # Ensure the temporary directory is removed upon exit
        trap "rm -rf $tmpdir" EXIT

        # Copy input files to the temporary directory
        cp "{input.bin}" "$tmpdir/"
        cp "{input.media}" "$tmpdir/"

        # Navigate to the temporary directory
        cd "$tmpdir"

        # Log the start of the carving process
        echo "Begin carving GEM for bin {params.binID} ..." >> "{log}"

        # Execute the carve command
        carve -g {params.carveme_medium} \
              -v \
              --mediadb "$(basename "{input.media}")" \
              --fbc2 \
              --solver scip \
              -o "$(basename "{input.bin}" .faa).xml" \
              "$(basename "{input.bin}")" >> "{log}" 2>&1

        # Log the completion of the carving process
        echo "Done carving GEM for bin {params.binID}." >> "{log}"

        # Move the resulting XML to the designated output directory
        mv "$(basename "{input.bin}" .faa).xml" "{output}"
        """


######################################################################################################
# Rule: modelVis -> one PDF + stats file per run (reading all GEMs together)
######################################################################################################
rule modelVis:
    input:
        gather_all_models
    output:
        text = os.path.join(OUT_PATH, "GEMs", "GEMs.stats"),
        plot = os.path.join(OUT_PATH, "Stats", "modelVis.pdf")
    conda:
        os.path.join(ENV_DIR, "metaGEMmod.yaml")
    log:
        os.path.join(OUT_PATH, "logs", "modelVis.log")
    shell:
        """
        mkdir -p $(dirname {output.text})
        mkdir -p $(dirname {output.plot})
        mkdir -p $(dirname {log})

        tmpdir=$(mktemp -d /dev/shm/modelvis.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT

        echo "Begin copying models to $tmpdir" &> "{log}"

        # Copy each model to a subdirectory named after the sample
        for f in {input}; do
            # sample is presumably the directory one level above GEMs
            sample=$(basename "$(dirname "$(dirname "$f")")")
            mkdir -p "$tmpdir/$sample"
            cp "$f" "$tmpdir/$sample/$(basename "$f")"
        done

        echo "Begin reading models and generating stats file..." >> "{log}"
        cd "$tmpdir"
        stats_file="GEMs.stats"
        rm -f "$stats_file"

        # Find all XML files in all subdirectories
        while IFS= read -r -d '' model; do
            sample_name=$(basename "$(dirname "$model")")
            bin_id="$(basename "$model" .xml)"  # e.g. bin.9
            mets=$(grep '<species ' "$model" | awk -F'id="' '{{print $2}}' \
                   | cut -d'"' -f1 | sort -u | wc -l)
            rxns=$(grep -c '</reaction>' "$model")
            genes=$(grep 'fbc:geneProduct=' "$model" | grep -vic 'spontaneous')
            
            echo "$sample_name $bin_id $mets $rxns $genes" >> "$stats_file"
        done < <(find . -name "*.xml" -print0)

        # Move the stats file into final output
        mv "$stats_file" "{output.text}"

        python {SCRIPTS_DIR}/modelVis.py \
            --stats "{output.text}" \
            --output "{output.plot}" >> "{log}" 2>&1
        """

######################################################################################################
# Rule: ECvis -> one directory of EC info per sample (aggregating all bin GEMs)
######################################################################################################
rule ECvis:
    input:
        gems = gatherbinmodels
    output:
        ec_dir = directory(os.path.join(OUT_PATH, "{sample}", "ecfiles"))
    log:
        os.path.join(OUT_PATH, "{sample}", "logs", "ECvis.log")
    shell:
        """
        # Create necessary directories
        mkdir -p {output.ec_dir}
        mkdir -p $(dirname "{log}")

        # Use /dev/shm for temporary storage
        tmpdir=$(mktemp -d /dev/shm/ECvis.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT

        echo "Extracting EC numbers from all GEMs for sample {wildcards.sample}..." &> {log}
        cp {input.gems} "$tmpdir/"
        cd "$tmpdir"

        mkdir ecfiles
        for model in *.xml; do
            grep -oP 'https://identifiers\\.org/ec-code/\\K[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+' "$model" \
                | sort | uniq -c > "ecfiles/${{model}}.ec"

            count=$(wc -l < "ecfiles/${{model}}.ec")
            echo "Model $model => $count unique EC numbers" >> {log}
        done

        mv ecfiles/* "{output.ec_dir}"/
        cd "{output.ec_dir}"
        awk '{{print $NF}}' *.ec | sort | uniq -c > EC.summary
        echo "EC.summary compiled for sample {wildcards.sample}" >> {log}
        """

######################################################################################################
# Rule: smetana -> single run per sample (all GEMs for that sample as input)
######################################################################################################

rule smetana:
    input:
        gatherbinmodels
    output:
        os.path.join(OUT_PATH, "{sample}", "SMETANA", "{sample}_detailed.tsv")
    params:
        media_db       = os.path.join(DB_DIR, "media_db.tsv"),
        smetana_media  = config["smetanaMedia"],
        smetana_solver = "cplex",
        smetana_flavor = "fbc2"
    conda:
        os.path.join(ENV_DIR, "smetana.yaml")
    log:
        os.path.join(OUT_PATH, "{sample}", "logs", "smetana.log")
    shell:
        """
        set -euo pipefail

        # Create necessary directories for output and logs
        mkdir -p "$(dirname "{output}")"
        mkdir -p "$(dirname "{log}")"

        # Create a unique temporary directory in /dev/shm with a recognizable prefix
        tmpdir=$(mktemp -d -p /dev/shm smetana_{wildcards.sample}_XXXXXX)
        echo "Temporary working directory: $tmpdir" >> "{log}"

        # Ensure the temporary directory is removed upon exit, regardless of success or failure
        trap "rm -rf $tmpdir" EXIT

        # Copy necessary input files to the temporary directory
        cp "{params.media_db}" "$tmpdir/"
        cp {input} "$tmpdir/"

        # Navigate to the temporary directory
        cd "$tmpdir"

        # (Optional) Load necessary modules if required
        module load math/CPLEX/22.11-GCCcore-10.2.0-Python-3.8.6 || true
        module load compiler/GCC/10.2.0 || true
        export LD_LIBRARY_PATH=${{CONDA_PREFIX}}/lib:$LD_LIBRARY_PATH

        # Log the start of the SMETANA process
        echo "Running SMETANA for sample {wildcards.sample} ..." >> "{log}"

        # Execute the SMETANA command with appropriate parameters
        smetana -o "{wildcards.sample}" \
                --flavor "{params.smetana_flavor}" \
                --mediadb "$(basename "{params.media_db}")" \
                -m "{params.smetana_media}" \
                --detailed \
                --solver "{params.smetana_solver}" \
                -v *.xml >> "{log}" 2>&1

        # Verify that the expected output file was created
        if [ ! -f "{wildcards.sample}_detailed.tsv" ]; then
            echo "Error: Expected output file {wildcards.sample}_detailed.tsv not found." >> "{log}"
            exit 1
        fi

        # Move the resulting TSV to the designated output directory
        mv "{wildcards.sample}_detailed.tsv" "{output}"

        # Log the completion of the SMETANA process
        echo "SMETANA done for sample {wildcards.sample}." >> "{log}"
        """

######################################################################################################
# Rule: memote -> one summary file per sample (covering all bin GEMs)
######################################################################################################
rule memote:
    input:
        gatherbinmodels
    output:
        out_dir = directory(os.path.join(OUT_PATH, "{sample}", "memote")),
        html = os.path.join(OUT_PATH, "{sample}", "memote", "{sample}.html")
    conda:
        os.path.join(ENV_DIR, "metaGEMmod.yaml")
    log:
        os.path.join(OUT_PATH, "{sample}", "logs", "memote.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})

        tmpdir=$(mktemp -d -p /dev/shm emote_{wildcards.sample}_XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT
        echo "Running memote for sample {wildcards.sample} ..." &> {log}

        cp {input} "$tmpdir/"
        cd "$tmpdir"

        # Compare all sample models
        echo "Running memote for all bin models ..." &>> {log}
        memote report diff \
            --filename {wildcards.sample}.html \
            {input} \
            --skip test_find_metabolites_produced_with_closed_bounds \
            --skip test_find_metabolites_consumed_with_closed_bounds \
            --skip test_find_metabolites_not_produced_with_open_bounds \
            --skip test_find_metabolites_not_consumed_with_open_bounds \
            --skip test_find_incorrect_thermodynamic_reversibility \
            >> "{log}" 2>&1
        
        # Move outputs to final out path
        echo "Moving outputs ..." &>> {log}
        # mv per_bin {wildcards.sample}.html {output.out_dir}/
        mv * {output.out_dir}/

        echo "Done." &>> {log}
        """
# rule GEM_Enzyme_Cluster:
#     input:
#         get_ec_dir, # Input directory with EC files
#         ec_summary = os.path.join(OUT_PATH, "ecfiles", "EC.summary"), # EC summary file
#     output:
#         enzymecluster = os.path.join(OUT_PATH, "enzymecluster", "enzymecluster.html")
#     conda:
#         os.path.join(ENV_DIR, "env.yaml")
#     log:
#         os.path.join(OUT_PATH, "logs", "enzymecluster.log")
#     shell:
#         """
#         python {SCRIPTS_DIR}/enzymecluster.py --ec_summary "{input.ec_summary}" \
#                                              --ec_dir "{input.get_ec_dir}" \
#                                              --output "{output.enzymecluster}" &> {log}
#         """
rule CompareBinsAndGems:
    input:
        bin_stats = glob.glob(directory(os.path.join(SAMPLES_PATH))), # Directory with bin stats files
        bins = glob.glob(directory(os.path.join(SAMPLES_PATH, "out"))), # Directory with protein FASTA files
        GEM_stats = "/mnt/lscratch/users/rparise/Thesis/Files/All_Samples/out/GEMs/GEMs.stats" # GEM stats file
    output:
        tsv = directory(os.path.join(OUT_PATH, "Stats/GEMs_vs_BinCompleteness"))
    conda:
        os.path.join(ENV_DIR, "env.yaml")
    log:
        os.path.join(OUT_PATH, "logs", "GEMs_vs_BinCompleteness.log")
    shell:
        """
        python {SCRIPTS_DIR}/ProteinsGemsVsBinCompleteness.py    --bin_stats "{input.bin_stats}" \
                                    --bins "{input.bins}" \
                                    --GEM_stats "{input.GEM_stats}" \
                                    --output "{output.tsv}" > {log} 2>&1
        """