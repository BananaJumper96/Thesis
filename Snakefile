outpath = "/mnt/lscratch/users/rparise/Thesis/Files/M11-07/Test_output"
scriptsDir = "/mnt/lscratch/users/rparise/Thesis/Scripts"
envDir = "/mnt/lscratch/users/rparise/Thesis/Envs"
databasesDir = "/mnt/lscratch/users/rparise/Thesis/Databases"

# Gather all protein fastas for each bin
def gatherbinmodels(wildcards):
    # Wait for the checkpoint to finish
    checkpoint_output = checkpoints.MakeBinProteinFastas.get(**wildcards).output.protbinfastaDir
    # Now the files should exist
    binIDs = glob_wildcards(os.path.join(checkpoint_output, "{binID}.faa")).binID
    # Return the paths to the GEM files
    return expand(os.path.join(outpath, "GEMs", "{binID}.xml"), binID=binIDs)
def get_detailed_tsv(wildcards):
    # Logic to generate binIDs
    checkpoint_output = checkpoints.MakeBinProteinFastas.get(**wildcards).output.protbinfastaDir
    binIDs = glob_wildcards(os.path.join(checkpoint_output, "{binID}.faa")).binID
    return expand(os.path.join(outpath, "SMETANA", "{binID}_detailed.tsv"), binID=binIDs)
def get_memote_dir(wildcards):
    # Logic to generate binIDs
    checkpoint_output = checkpoints.MakeBinProteinFastas.get(**wildcards).output.protbinfastaDir
    binIDs = glob_wildcards(os.path.join(checkpoint_output, "{binID}.faa")).binID
    return expand(os.path.join(outpath, "memote", "{binID}"), binID=binIDs)

rule all:
    input:
        gatherbinmodels,
        get_memote_dir, 
        get_detailed_tsv, # Smetana rule
        ec_dir = directory(os.path.join(outpath, "ecfiles")),
        text = os.path.join(outpath, "GEMs", "GEMs.stats"),
        plot = os.path.join(outpath, "Stats", "modelVis.pdf"),
        
        
checkpoint MakeBinProteinFastas:
    input:
        cont2bin = "/mnt/lscratch/users/rparise/Thesis/Files/M11-07/bins.contigs.gz",
        cont2orf = "/mnt/lscratch/users/rparise/Thesis/Files/M11-07/annotation.filt.contig2ID.tsv.gz",
        protein_fasta = "/mnt/lscratch/users/rparise/Thesis/Files/M11-07/annotation.faa.gz"
    output:
        protbinfastaDir = directory(os.path.join(outpath, "protein_bins"))
    log:
        os.path.join(outpath, "logs/MakeBinProteinFastas.log")
    conda:
        os.path.join(envDir, "env.yaml")
    shell:
        """
        python {scriptsDir}/getbinprotfastas.py --cont2bin {input.cont2bin} \
                                                --cont2orf {input.cont2orf} \
                                                --protein_fasta {input.protein_fasta} \
                                                --output_dir {output.protbinfastaDir} &> {log}
        """
# Rule to generate GEMs from protein FASTA files using CarveMe
rule carveme:
    input:
        bin = os.path.join(outpath,"protein_bins/{binID}.faa"),
        media = os.path.join(databasesDir, "media_db.tsv")   # Media database file for CarveMe
    output:
        os.path.join(outpath,"GEMs/{binID}.xml")  # Output GEM in SBML format
    params:
        binID = lambda wildcards: wildcards.binID  # Extract bin ID from the input filename
    threads: 4
    message:
        """
        Running CarveMe
        """
    conda:
        os.path.join(envDir, "env.yaml")
    log:
        os.path.join(outpath, "logs/carveme_{binID}.log")
    shell:
        """
        set -euo pipefail

        # Create a temporary directory
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        # Copy input files to the temporary directory
        cp "{input.bin}" "$tmpdir/"
        cp "{input.media}" "$tmpdir/"

        cd "$tmpdir"

        echo "Begin carving GEM ..." >> "{log}"

        # Run CarveMe to generate the GEM
        carve -g "M8" \
            -v \
            --mediadb "$(basename "{input.media}")" \
            --fbc2 \
            -o "$(basename "{input.bin}" .faa).xml" \
            "$(basename "{input.bin}")" >> "{log}" 2>&1

        echo "Done carving GEM." >> "{log}"

        # Ensure the output directory exists
        mkdir -p "$(dirname "{output}")"

        # Move the generated GEM to the output directory
        if [ -f "$(basename "{input.bin}" .faa).xml" ]; then
            mv "$(basename "{input.bin}" .faa).xml" "{output}"
        else
            echo "Error: GEM file not generated." >> "{log}"
            exit 1
        fi
        """
# Rule to generate statistics and visualize GEMs
rule modelVis:
    input: 
        gatherbinmodels
    output: 
        text = os.path.join(outpath, "GEMs", "GEMs.stats"),
        plot = os.path.join(outpath, "Stats", "modelVis.pdf")
    message:
        """
        Generate bar plot with GEMs generated across samples and density plots showing number of 
        unique metabolites, reactions, and genes across GEMs.
        """
    conda:
        os.path.join(envDir, "env.yaml")
    log:
        os.path.join(outpath, "logs/modelVis.log")
    shell:
        """
        # Copy the GEMs to the scratch directory
        mkdir -p /dev/shm/
        cp {input} /dev/shm/
        cd /dev/shm/
        ls &> {log}
        echo -e "\nBegin reading models ... \n" &>> {log}
        
        # Initialize the statistics file
        stats_file="GEMs.stats"
        rm -rf "$stats_file" &>> {log}
        echo starting Loop &>> {log}

        # Loop over each GEM file and extract statistics
        for model in *.xml; do 
            echo "Processing model: $model ... " &>> {log}
            id="${{model%.xml}}"  # Extract model ID by removing .xml extension
            echo "Model ID: $id" &>> {log}
            mets=$(grep '<species ' "$model" 2>> {log} | awk -F'id="' '{{print $2}}' 2>> {log} | cut -d'"' -f1 2>> {log} | sort -u 2>> {log} | wc -l 2>> {log})
            echo "Model: $id has $mets mets ... " &>> {log}
            rxns=$(grep -c "</reaction>" "$model" 2>> {log}) 
            echo "Model: $id has $rxns reactions ... " &>> {log}
            genes=$(grep "fbc:geneProduct=" "$model" 2>> {log} | grep -vic "spontaneous" 2>> {log})
            echo "Model: $id has $genes genes ... " &>> {log}
            echo "Model: $id has $mets mets, $rxns reactions, and $genes genes ... " >> {log}
            echo "$id $mets $rxns $genes" >> "$stats_file" 2>> {log} # Append statistics to file 
        done

        echo -e "\nDone generating GEMs.stats summary file, moving to output folder and running modelVis.py script ... " &>> {log}

        # Move the statistics file to the output directory
        mv "$stats_file" "{output.text}"
        echo "Stats file moved to output directory." &>> {log}
        echo $stats_file &>> {log}
        # Run Python script to generate visualizations
        python {scriptsDir}/modelVis.py --stats "{output.text}" \
                                        --output "{output.plot}" &>> {log}
        
        echo "Done." &>> {log}
        """

# Rule to extract EC numbers from GEMs and create summary statistics
rule ECvis:
    input: 
        gems = gatherbinmodels  # Input list of GEM model files
    output:
        ec_dir = directory(os.path.join(outpath, "ecfiles"))  # Output directory for EC files
    log:
        os.path.join(outpath, "logs/ECvis.log")  # Log file
    message:
        """
        Extract EC information from GEMs.
        """
    shell:
        """
        set -euo pipefail  # Enable strict error handling

        # Create a unique temporary directory in /dev/shm
        tmpdir=$(mktemp -d /dev/shm/ec_extraction.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT  # Ensure the temporary directory is removed on exit

        echo -e "\\nBegin reading models ... \\n" &>> {log}

        # Copy the GEMs to the temporary directory
        cp {input.gems:q} "$tmpdir/"

        cd "$tmpdir"

        # Create directory for EC files
        mkdir -p ecfiles

        for model in *.xml; do
            # Extract EC numbers from the GEM and write to a unique file
            echo "Reading EC numbers in model $model ..." &>> {log}
            grep -oP 'https://identifiers\\.org/ec-code/\\K[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+' "$model" | \
                sort | uniq -c > "ecfiles/${{model}}.ec"

            # Count unique EC numbers
            ECNUM=$(wc -l < "ecfiles/${{model}}.ec")

            echo "Model $model has $ECNUM unique EC numbers." &>> {log}
        done

        echo -e "\\nMoving ecfiles folder back to {output.ec_dir}" &>> {log}
        mv "ecfiles" "{output.ec_dir}"

        cd "{output.ec_dir}"

        echo -e "\\nCreating sorted unique file EC.summary for easy EC inspection ..." &>> {log}

        # Combine all EC files into a summary
        cat *.ec | awk '{{print $NF}}' | sort | uniq -c > EC.summary

        # Display the EC summary
        cat EC.summary &>> {log}

        echo "Done." &>> {log}
        """
# Rule to generate interaction visualization from SMETANA results
# rule interactionVis:
#     input:
#         smetana_results = expand(os.path.join(outpath, "SMETANA", "{binID}_detailed.tsv"), binID=binIDs)
#     output:
#         stats = os.path.join(outpath, "Stats", "sampleMedia.stats")
#     params:
#         binID = lambda wildcards: wildcards.binID  # Extract bin ID from the input filename
#     log:
#         os.path.join(outpath, "logs", "interactionVis.log")
#     shell:
#         """
#         set -euo pipefail  # Enable strict error handling

#         # Change to the SMETANA results directory
#         cd "{outpath}/SMETANA"

#         # Combine all SMETANA result files, excluding "community" lines
#         grep -h -v "community" *.tsv > smetana.all

#         # Extract unique media types
#         cut -f2 smetana.all | sort | uniq > media.txt

#         # List sample IDs from the input files
#         ls *.tsv | sed "s/_detailed.tsv$//" > samples.txt

#         # Generate statistics of media usage per sample
#         {{
#             while read -r sample; do 
#                 echo -n "$sample "
#                 while read -r media; do 
#                     var=$(grep -F "$sample" smetana.all | grep -F -c "$media")
#                     echo -n "$var " 
#                 done < media.txt
#                 echo ""
#             done < samples.txt
#         }} > "{output.stats}"

#         echo "Interaction visualization statistics generated at {output.stats}" &>> {log}
#         """
# Rule to run SMETANA for metabolic interaction analysis per sample
rule smetana:
    input:
        gatherbinmodels
    output:
        detailed_tsv = os.path.join(outpath, "SMETANA", "{binID}_detailed.tsv")
    params:
        smetana_media = "M8",  # Replace with your media
        smetana_solver = "cplex",    # Replace with your solver
        smetana_flavor = "fbc2", # Replace with your model flavor
        media_db = os.path.join(databasesDir, "media_db.tsv"), # Media database file for SMETANA
        binID = lambda wildcards: wildcards.binID  # Extract bin ID from the input filename
    conda:
        os.path.join(envDir, "env_gurobi_02.yaml") # Use a separate environment with Gurobi
    log:
        os.path.join(outpath, "logs", "{binID}_smetana.log") # Log file
    shell:
        """
        module load math/CPLEX/22.11-GCCcore-10.2.0-Python-3.8.6
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        mkdir -p "$(dirname "{output.detailed_tsv}")"

        echo -e "\\nRunning SMETANA for sample {wildcards.binID} ... " >> "{log}"

        cp "{params.media_db}" "$tmpdir/" >> "{log}" 2>&1
        cp {input:q} "$tmpdir/" >> "{log}" 2>&1

        cd "$tmpdir"
        smetana -o "{wildcards.binID}" \
            --flavor "{params.smetana_flavor}" \
            --mediadb "$(basename "{params.media_db}")" \
            -m "{params.smetana_media}" \
            --detailed \
            --solver "{params.smetana_solver}" \
            -v *.xml >> "{log}" 2>&1

        cp "{wildcards.binID}_detailed.tsv" "{output.detailed_tsv}" >> "{log}" 2>&1

        echo "SMETANA analysis for sample {wildcards.binID} completed." >> "{log}"
        """
# Rule to run Memote for quality control of GEMs
rule memote:
    input:
        gem_file = os.path.join(outpath, "GEMs", "{binID}.xml")
    output:
        memote_dir = directory(os.path.join(outpath, "memote", "{binID}"))
    params:
        skip_tests = [
            "--skip", "test_find_metabolites_produced_with_closed_bounds",
            "--skip", "test_find_metabolites_consumed_with_closed_bounds",
            "--skip", "test_find_metabolites_not_produced_with_open_bounds",
            "--skip", "test_find_metabolites_not_consumed_with_open_bounds",
            "--skip", "test_find_incorrect_thermodynamic_reversibility"
        ],
        binID = lambda wildcards: wildcards.binID  # Extract bin ID from the input filename
    conda:
        os.path.join(envDir, "env.yaml")
    log:
        os.path.join(outpath, "logs", "{binID}_memote.log")
    shell:
        """
        set -euo pipefail  # Enable strict error handling

        tmpdir=$(mktemp -d /dev/shm/memote_{wildcards.binID}_XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT

        echo -e "\\nRunning Memote for GEM {wildcards.binID} ... " &>> {log}

        cp "{input.gem_file}" "$tmpdir/"
        cd "$tmpdir"

        memote report snapshot \
            {params.skip_tests} \
            --filename "{wildcards.binID}.html" \
            "{input.gem_file}" &>> {log}

        memote run \
            {params.skip_tests} \
            "{input.gem_file}" &>> {log}

        mv "result.json.gz" "{wildcards.binID}.json.gz"

        mkdir -p "{output.memote_dir}"
        mv *.gz *.html "{output.memote_dir}/"

        echo "Memote analysis for GEM {wildcards.binID} completed." &>> {log}
        """
rule gathermodels:
    input:
        gatherbinmodels
    output:
        os.path.join(outpath, "gathermodels.txt")
    shell:
        """
        touch {output}
        """
